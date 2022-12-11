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
rs_train <- lapply(1:10, function(i) {temp <- rs_resample[[i]]
                        list(brain_train=temp$pca_train$brain_train_reduced,
                        cbcl_train=temp$cbcl_train)})                  
saveRDS(rs_train, "rs_train.rds")

### extract the brain and CBCL features in the 10 test sets
rs_test <- lapply(1:10, function(i) {temp <- rs_resample[[i]]
                        # the brain PCs in the test sets are calculated by applying the eigenvectors in the training sets
                        brain_test <- as.matrix(temp$brain_test) %*% temp$pca_train$rotation[, 1:100]
                        list(brain_test=brain_test,
                             cbcl_test=temp$cbcl_test)})
saveRDS(rs_test,"rs_test.rds")

### extract the eigenvectors of the brain PCs in the 10 training sets
rs_rotation <- lapply(1:10, function(i) {temp <- rs_resample[[i]]
                        brain_rotation <- temp$pca_train$rotation[, 1:100]})
saveRDS(rs_rotation,"rs_rotation.rds")

rm(rs_resample)


#########################################################
############## 2. grid search of penalty parameters #####
#########################################################

# grid search of penalty parameters in ABCD training set

grid.abcd <- lapply(1:10, function(i) {
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

########################################################
############## 3. Fit the sCCA model  ##################
########################################################


### fit the CCA model in 10 train-test splits and Generation R

rs_train_test_abcd <- lapply(1:10, function(i) {
  
  brain_train <- rs_train[[i]]$brain_train
  cbcl_train <- rs_train[[i]]$cbcl_train
  brain_test <- rs_test[[i]]$brain_test
  cbcl_test <- rs_test[[i]]$cbcl_test
  brain_pen <- grid.abcd[[i]][[1]][[1]]
  cbcl_pen <- grid.abcd[[i]][[1]][[2]]
  
  # fit the SCCA model in the training set
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  # permutation test of the canonical correlations in the training set
  perm_abcd_train <- permutation_test(cbcl_train, brain_train, nperm=999, cbcl_pen, brain_pen,8,res.abcd$cors)
  
  # project the CCA weights in the ABCD test set
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd,8)
  # permutation test of the canonical correlations in the test set
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=1999,res.abcd,abcd.test)
  
  # project the CCA weights in Generation R dara
  cor.abcdTogenr <- test_project_weights(brain_genr, cbcl_genr, res.abcd, 8)
  
  perm_abcdTogenr <- permutation_test_testset(cbcl_genr,brain_genr,nperm=1999,res.abcd,abs(cor.abcdTogenr))
  
  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr$pval.perm))
})










