#####################################################################
############## ABCD: grid search, sparse CCA, bootstrapping #########
#####################################################################


########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'fmsb', 'NbClust','stats','dendextend','cluster',
              'fpc', 'e1071', 'plotly', 'MASS', 'ggradar2','tidyr','introdatavizyes','Hmisc',
              'pheatmap','RColorBrewer','circlize')


lapply(packages, require, character.only = TRUE)

########################################################
############## 1. Read the data ########################
########################################################

### read the data
rs_resample <- readRDS("train_test_split.rds")

### extract the brain and CBCL features in the 10 training sets
rs_train <- lapply(1:30, function(i) {temp <- rs_resample[[i]]
                        list(brain_train=temp$pca_train$brain_train_reduced,
                        cbcl_train=temp$cbcl_train)})                  
saveRDS(rs_train, "rs_train.rds")

### extract the brain and CBCL features in the 10 test sets
rs_test <- lapply(1:30, function(i) {temp <- rs_resample[[i]]
                        # the brain PCs in the test sets are calculated by applying the eigenvectors in the training sets
                        brain_test <- as.matrix(temp$brain_test) %*% temp$pca_train$rotation[, 1:100]
                        list(brain_test=brain_test,
                             cbcl_test=temp$cbcl_test)})
saveRDS(rs_test,"rs_test.rds")

### extract the eigenvectors of the brain PCs in the 10 training sets
rs_rotation <- lapply(1:30, function(i) {temp <- rs_resample[[i]]
                        brain_rotation <- temp$pca_train$rotation[, 1:100]})
saveRDS(rs_rotation,"rs_rotation.rds")

rm(rs_resample)


#########################################################
############## 2. grid search of penalty parameters #####
#########################################################

# grid search of penalty parameters in ABCD training set

grid.abcd <- lapply(1:30, function(i) {
  rsfmri_resample <- rs_train[[i]]
  brain_train <- rsfmri_resample$brain_train
  cbcl_train <- rsfmri_resample$cbcl_train
  
  # split the ABCD training set into training (80%) and validation (20%) sets 100 times
  cv.resample <- CV_sampling(cbcl_train, brain_train, 100, 0.8)
  
  x_pen <- seq(0.1, 1.0, length.out = 10) # brain
  y_pen <- seq(0.1, 1.0, length.out = 10) # cbcl
  
  grid.abcd <- grid.search.cor.Testset(cv.resample$brain_train, cv.resample$cbcl_train,
                                       cv.resample$brain_test, cv.resample$cbcl_test,
                                       x_pen,y_pen,nsample=100)
})

saveRDS(grid.abcd,"grid.abcd.rds")

###### you can check the loadings in your cases and decide whether you want to constrain the penalty parameters
###### the CBCL loadings were constrained to be larger than 0.5 in the current study to improve the intepretability 
grid.abcd.fmriprep.constrained <- lapply(seq_along(join.abcd), function(i) {
  max_cor <- max(join_grid[[i]][[2]][, 5:10])
  idx <- which(join_grid[[i]][[2]] == max_cor, arr.ind = T)
  penalty_brain <- idx[1,1] * 0.1
  penalty_cbcl <- idx[1,2] * 0.1
  list(penalty_brain=penalty_brain, penalty_cbcl=penalty_cbcl)
})


########################################################
############## 3. Fit the sCCA model  ##################
########################################################


### fit the CCA model in 30 train-test splits and Generation R

rs_train_test_abcd <- lapply(1:30, function(i) {
  
  brain_train <- rs_train[[i]]$brain_train
  cbcl_train <- rs_train[[i]]$cbcl_train
  brain_test <- rs_test[[i]]$brain_test
  cbcl_test <- rs_test[[i]]$cbcl_test
  brain_pen <- grid.abcd.fmriprep.constrained[[i]]$penalty_brain
  cbcl_pen <- grid.abcd.fmriprep.constrained[[i]]$penalty_cbcl
  
  # fit the SCCA model in the training set
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  # permutation test of the canonical correlations in the training set
  perm_abcd_train <- permutation_test(cbcl_train, brain_train, nperm=1999, cbcl_pen, brain_pen,3,res.abcd$cors) # based on the scree plot of covariance explained
  
  # project the CCA weights in the ABCD test set
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd,8)
  # permutation test of the canonical correlations in the test set
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=1999,res.abcd,abcd.test[1:3])
  
  # project the CCA weights in Generation R dara
  cor.abcdTogenr <- test_project_weights(brain_genr, cbcl_genr, res.abcd, 8)
  
  perm_abcdTogenr <- permutation_test_testset(cbcl_genr,brain_genr,nperm=1999,res.abcd,abs(cor.abcdTogenr))
  
  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr$pval.perm))
})


##################################################
############## multiple testing correction #######
##################################################
lapply(1:30, function(x) {p.adjust(rs_train_test_abcd[[x]]$abcd.train.perm, "fdr")})
lapply(1:30, function(x) {p.adjust(rs_train_test_abcd[[x]]$abcd.test.perm, "fdr")})
lapply(1:30, function(x) {p.adjust(rs_train_test_abcd[[x]]$genr.perm, "fdr")})


