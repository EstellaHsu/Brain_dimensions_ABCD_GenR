#################################################################
############## ABCD: training test with residualization #########
#################################################################


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
setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/retraintest")
### read the data
rs_resample <- lapply(c(1:3,5:11), function(i) { path <- paste0("train_test_split_", i, ".rds")
                           rs_sample <- readRDS(path)})

rm(rs_resample)

rs_train <- lapply(1:10, function(i) {temp <- rs_resample[[i]]
                        list(brain_train=temp$pca_train$brain_train_reduced,
                        cbcl_train=temp$cbcl_train)})
saveRDS(rs_train, "rs_train.rds")
rs_test <- lapply(1:10, function(i) {temp <- rs_resample[[i]]
                        brain_test <- as.matrix(temp$brain_test) %*% temp$pca_train$rotation[, 1:100]
                        list(brain_test=brain_test,
                             cbcl_test=temp$cbcl_test)})
saveRDS(rs_test,"rs_test.rds")
rs_rotation <- lapply(1:10, function(i) {temp <- rs_resample[[i]]
                        brain_rotation <- temp$pca_train$rotation[, 1:100]})
saveRDS(rs_rotation,"rs_rotation.rds")

rm(rs_resample)

p1 <- readRDS("train_test_split_1.rds")
f1 <- p1$pca_train$brain_train_reduced
f1_cbcl <- p1$cbcl_train


lapply(rs_test, function(x) {nrow(x$cbcl_test)})
################################################
############## 2. grid search ##################
################################################
####### Step 1: train penalty parameters 
##### weighted PCA
grid.abcd.wholesample <- lapply(1:10, function(i) {
  rsfmri_resample <- rs_train[[i]]
  brain_train <- rsfmri_resample$brain_train
  cbcl_train <- rsfmri_resample$cbcl_train
  cv.resample <- CV_sampling(cbcl_train, brain_train, 100, 0.8)
  
  x_pen <- seq(0.1, 1.0, length.out = 10) # brain
  y_pen <- seq(0.1, 1.0, length.out = 10)
  
  grid.abcd <- grid.search.cor.Testset(cv.resample$brain_train, cv.resample$cbcl_train,
                                       cv.resample$brain_test, cv.resample$cbcl_test,
                                       x_pen,y_pen,nsample=100)
})

saveRDS(grid.abcd.wholesample,"grid.abcd.wholesample.rds")

# 0.5, 0.8
# 0.5, 0.5
# 0.5, 0.3/0.2
# 0.5, 0.5
# 0.5, 0.6
# 0.5, 0.3
# 0.5, 0.5
# 0.5, 0.5
# 0.7/0.5, 0.2
# 0,5, 0.3

########################################################
############## 3. Fit the sCCA model  ##################
########################################################
res.abcd.5 <- CCA(x=rs_train[[5]]$brain_train, z=rs_train[[5]]$cbcl_train, penaltyx = 0.6, penaltyz = 0.5, typex="standard", typez="standard",
                          niter = 20, K=8)
res.abcd.5

# visualize the loadings:
cbcl_loading <- res.abcd.1$v

#cbcl_loading[, 3:5] <- -cbcl_loading[, 3:5] 
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")
corrplot(t(cbcl_loading)[1:3,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 1.6, number.cex=1)

### permutation test
perm_abcd1 <- permutation_test(f1_cbcl, f1, nperm=999, 0.2, 0.5, 8, res.abcd.siemens.1$cors)

df.cor.2 <- data.frame(cor = perm_abcdtest10$cor.perm[, 2])
df.cor.1 <- data.frame(cor = perm_abcdtest10$cor.perm[, 1])
df.cor.3 <- data.frame(cor = perm_abcdtest10$cor.perm[, 3])

pdf("permutation_test_cv3.pdf", width=5, height=3)
ggplot(df.cor.3, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "#afe0a9", color = "black", alpha = 0.6)+
  geom_vline(xintercept=0.13, linetype="dashed", color = "red", size = 1) +
  ggtitle("Permutation correlations of the 1st canonical variate") + 
  labs(x = "canonical correlations") + 
  #annotate(geom="text", x=0.0, y=200, label="r = 0.13", size = , color="black") + 
  xlim(c(-0.2, 0.2)) + annotate(geom="text", x=0.0, y=180, label="P(FDR) < 0.01", 
                                fontface = 'italic', size = 5, color="black") + 
  theme_bw()
dev.off()

cor2 <- ggplot(df.cor.2, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "lightblue", color = "black", alpha = 0.8)+
  geom_vline(xintercept=0.07, linetype="dashed", color = "red", size = 1) + 
  annotate(geom="text", x=0.13, y=250, label="r = 0.16", size = 5, color="black") + xlim(c(-0.2, 0.2))+
  annotate(geom="text", x=0.13, y=225, label="p(fdr) < 0.001", fontface = 'italic', size = 5, color="black") +
  ggtitle("Permutation correlations of the 2nd canonical variate") + 
  labs(x = "canonical correlations") + 
  theme_bw()
cor3 <- ggplot(df.cor.3, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "limegreen", color = "black", alpha = 0.6)+
  geom_vline(xintercept=0.09, linetype="dashed", color = "red", size = 1) + 
  annotate(geom="text", x=0.13, y=250, label="r = 0.17", size = 5, color="black") + xlim(c(-0.2, 0.2))+
  annotate(geom="text", x=0.13, y=225, label="p(fdr) < 0.001", fontface = 'italic', size = 5, color="black") +
  ggtitle("Permutation correlations of the 3rd canonical variate") + 
  labs(x = "canonical correlations") + 
  theme_bw()

grid.arrange(cor1, cor2, cor3, ncol =3)  

########################################################
############## Covariance explained  ###################
########################################################
vardf <- VarianceExplain(f7_syn, f7_syn_cbcl, res.wei.syn.abcd7, 8) 
colnames(vardf) <- c("principal_components", "Covariance_explained")
ggplot(vardf, aes(x=principal_components, y=Covariance_explained))+
  geom_point(size = 3) + theme_bw()

############################################################
############## Replication within ABCD  ####################
############################################################
test_abcd <- as.matrix(p1$brain_test) %*% p1$pca_train$rotation[, 1:100]
dim(test_abcd)
# direct mapping
cor.abcdtest <- test_project_weights(test_abcd, p1$cbcl_test,res.abcd.1, 8)
cor.abcdtest
perm_abcdtest <- permutation_test_testset(p1$cbcl_test, test_abcd,nperm=999, 
                                                  res.abcd.1,cor.abcdtest)
p.adjust(perm_abcdtest1$pval.perm, "fdr")


###############################################################
############## Replication in Generation R  ####################
###############################################################
setwd("/Users/estella/Desktop/ABCD_download/data/Generation_R")
feature_brain_centered <- readRDS("feature_brain_centered_genr.rds")
test.pca.genr <- as.matrix(feature_brain_centered) %*% as.matrix(p1$pca_train$rotation[, 1:100])

# direct mapping
cor.abcdTogenr <- test_project_weights(test.pca.genr, cbcl_syn_genr, res.abcd.1, 6)
cor.abcdTogenr

########################################
############ automatic process
########################################
# 0.5, 0.8 1
# 0.5, 0.5 2
# 0.5, 0.3/0.2 3
# 0.5, 0.5 4
# 0.5, 0.6 5
# 0.5, 0.3 6
# 0.5, 0.5 7
# 0.5, 0.5 8
# 0.7/0.5, 0.2 9
# 0,5, 0.3 10
cbcl_pen <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.7,0.5)
brain_pen <- c(0.8,0.5,0.3,0.5,0.6,0.3,0.5,0.5,0.2,0.3)

setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/retraintest")
rs_train <- readRDS("rs_train.rds")
rs_test <- readRDS("rs_test.rds")
rs_rotation <- readRDS("rs_rotation.rds")

rs_train_test_abcd <- lapply(1:10, function(i) {
  
  brain_train <- rs_train[[i]]$brain_train
  cbcl_train <- rs_train[[i]]$cbcl_train
  brain_test <- rs_test[[i]]$brain_test
  cbcl_test <- rs_test[[i]]$cbcl_test
  brain_pen <- brain_pen[i]
  cbcl_pen <- cbcl_pen[i]
  
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  perm_abcd_train <- permutation_test(cbcl_train, brain_train, nperm=999, cbcl_pen, brain_pen, 8, res.abcd$cors)
  
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd, 8)
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=999,res.abcd,abcd.test)
  
  cor.abcdTogenr <- test_project_weights(brain_genr_str, cbcl_genr_str, res.abcd, 8)
  perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr_str,brain_genr_str,nperm=999, 
                                                    res.abcd,abs(cor.abcdTogenr))
  
  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr_total$pval.perm))
})


cbcl_loading <- abs(rs_train_test_abcd[[7]]$abcd.train$v)
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")

corrplot(t(cbcl_loading)[1:5,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 2, number.cex=1)

vardf <- VarianceExplain(rs_train[[5]]$brain_train, rs_train[[5]]$cbcl_train, res.abcd.5, 8) 
colnames(vardf) <- c("principal_components", "Covariance_explained")
ggplot(vardf, aes(x=principal_components, y=Covariance_explained))+
  geom_point(size = 3) + theme_bw()

saveRDS(rs_train_test_abcd,"rs_train_test_abcd.rds")



##############################################################################
############ use split 7 to further explore the canonical loadings ###########
##############################################################################

rsbrain_loading <- abs(rs_train_test_abcd[[7]]$abcd.train$u)
rsbrain_loading <- rsbrain_loading[, 1:3]
rownames(rsbrain_loading) <- paste0("PC",1:nrow(rsbrain_loading))
colnames(rsbrain_loading) <- rownames(c("CV1","CV2","CV3"))

rs_brain_loadings <- rsbrain_loading[rowSums(rsbrain_loading) != 0, ]


###########################################################################
###### bootstrap of the best split: stability check of cbcl and the brain
###########################################################################

# split 7 as an example
setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/retraintest")
### read the data
rs_example <- readRDS("train_test_split_7.rds")
brain_train <- rs_example$pca_train$brain_train_reduced
cbcl_train <- rs_example$cbcl_train
brain_test <- as.matrix(rs_example$brain_test) %*% rs_example$pca_train$rotation[, 1:100]
cbcl_test <- rs_example$cbcl_test 
res.abcd.7 <- CCA(x=brain_train, z=cbcl_train, penaltyx = 0.5, penaltyz = 0.5, typex="standard", typez="standard",
                  niter = 20, K=8)
res.abcd.7


######## bootstrap of the CVs and correlations
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
n_boot <- 1000
boot_abcd <- foreach(1:n_boot, .packages="PMA") %dopar% {
  # boot samples
  N <- nrow(cbcl_train)
  idx <- sample.int(N, N, replace = TRUE) 
  # index the boot samples
  brain_train_boot <- brain_train[idx, ]
  cbcl_train_boot <- cbcl_train[idx, ]
  # run cca
  abcd_boot <- CCA(x=brain_train_boot, z=cbcl_train_boot, penaltyx = 0.5, penaltyz = 0.5, 
                   typex="standard", typez="standard",niter = 20, K=8)
  # extract loadings
  cbcl_loading <- abcd_boot$v
  brain_loading <- abcd_boot$u
  # test on the test set
  cor.abcd.test_boot <- test_project_weights(brain_test, cbcl_test, abcd_boot, 8)
  
  list(cbcl_loading=cbcl_loading, brain_loading=brain_loading, 
       cors_train=abcd_boot$cors, cors_test=cor.abcd.test_boot)
}
stopCluster(cl)

###### reorder the CVs of cbcl and brain
load_std <- res.abcd.7$u
cv_reorder_boot_abcd  <- lapply(seq_along(boot_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:4) {
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:4, function(z) {cor(load_std[, i], boot_abcd[[j]]$brain_loading[, z])})
    if(max(abs(cor_cv)) < 0.5){ # sometimes the correlation is too low
      idx  <- i
    } else {
      idx <- which(abs(cor_cv) == max(abs(cor_cv)))
    }
    
    if(!idx %in% idx_v) {
      idx_v <- c(idx_v, idx) # sometimes it will overlap
    } else {
      idx_v <- c(idx_v, i) 
    }
  }
  
  idx_v[duplicated(idx_v)] <- c(1:4)[!c(1:4) %in% idx_v]
  list(idx_reorder=idx_v)
  
})



# training set
boot_cor_train <- lapply(boot_abcd, function(x) {x$cors_train[1:4]})
boot_cor_train <- do.call(rbind,boot_cor_train)
boot_cor_train_reorder <- lapply(1:n_boot, function(i) {boot_cor_train[i, cv_reorder_boot_abcd[[i]]$idx_reorder]})
boot_cor_train_reorder <- do.call(rbind, boot_cor_train_reorder)
colnames(boot_cor_train_reorder) <- paste0("CV",1:4)
meantrain_boot <- colMeans(boot_cor_train_reorder)
sdtrain_boot <- apply(boot_cor_train_reorder,2,sd)

# test set
boot_cor_test <- lapply(boot_abcd, function(x) {x$cors_test[1:4]})
boot_cor_test <- do.call(rbind,boot_cor_test)
boot_cor_test_reorder <- lapply(1:n_boot, function(i) {boot_cor_test[i, cv_reorder_boot_abcd[[i]]$idx_reorder]})
boot_cor_test_reorder <- do.call(rbind,boot_cor_test_reorder)
colnames(boot_cor_test_reorder) <- paste0("CV",1:4)
meantest_boot <- colMeans(boot_cor_test_reorder)
sdtest_boot <- apply(boot_cor_test, 2, sd)


fig_boot_syn <- data.frame(Train_Test = rep(c("Training set","Test set"),each=4), meancor = c(meantrain_boot, meantest_boot),
                           sdcor = c(sdtrain_boot, sdtest_boot), CV=rep(paste0("CV",1:4),2))

########### violin plot with error bars
#########################
#train
df1 <- as.data.frame(boot_cor_train_reorder[,1:3])
names(df1) <- paste0("CV", 1:3)
df1$traintest <- "Train"
#test
df2 <- as.data.frame(boot_cor_test_reorder[,1:3])
names(df2) <- paste0("CV", 1:3)
df2$traintest <- "Test"

# combine train and test
df3 <- rbind(df1,df2)
# reshape the data frame
df_long <- gather(df3, CV, cors, CV1:CV3, factor_key=TRUE)
df_long$cors <- abs(df_long$cors)


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


p3 <- ggplot(df_long, aes(x = CV, y = abs(cors), fill = traintest)) +
  geom_point(aes(col=traintest), position=position_jitter(height=.05, width=.25),alpha = .3)+ 
  introdataviz::geom_split_violin(alpha = .8, trim = FALSE) +
  stat_summary(fun.data = data_summary,geom="pointrange", show.legend = TRUE, 
               position = position_dodge(.3)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text = element_text(size=13))+
  scale_x_discrete(name = "Canonical Variates", labels = paste0("CV", 1:3)) +
  scale_y_continuous(name = "Canonical correlations") +
  scale_fill_manual(name="Training Test Set", values = c("#fbb172", "#d45365")) + 
  scale_color_manual(name="Training Test Set", values = c("#fbb172", "#d45365")) + 
  theme_bw()

p3
ggsave("train_test_boot_split7_3_violin_reorder.pdf", width = 8, height = 6)

boot_abcd[[31]]$cbcl_loading

######################################################
####### visualization of the CIs of cbcl loadings ####
######################################################
cbclnames <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule Breaking","Aggression")
# training set
n_boot <- 1000
boot_cbcl <- lapply(1:n_boot, function(i) {
  boot_abcd[[i]]$cbcl_loading[,cv_reorder_boot_abcd[[i]]$idx_reorder]})

cbcl_sign <- sign(res.abcd.7$v)


# cv1: attention, cv2: rule breaking, cv3:withdrawn, all +1

idx <- c(6,8,2) # this is the position of the most important cbcl syndrome scale (attention, aggression, withdrawn)
df_cbcl <- lapply(1:3, function(i) {
  cbcl_cv <- do.call(cbind, lapply(boot_cbcl, function(x){x[, i]}))
  cbcl_cv <- t(abs(cbcl_cv))
  cbcl_cv <- as.data.frame(cbcl_cv)
  # adjust the signs 
  for (j in 1:nrow(cbcl_cv)){
    if(sign(cbcl_cv[j, idx[i]]) != 1){
      cbcl_cv[j,] <- -cbcl_cv[j, ]
    }
  }
  names(cbcl_cv) <- cbclnames
  cbcl_cv$CV <- paste0("CV", i)
  return(cbcl_cv)})

df_cbcl <- do.call(rbind,df_cbcl)
df_long_cbcl <- gather(df_cbcl, cbclsyndromes, loadings, cbclnames, factor_key = TRUE)


############ combine three CVs together
# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND treatment :
pdf(file="cbcl_loadings_boot_newsign.pdf", width=13, height=5)
myplot <- boxplot(loadings ~ CV*cbclsyndromes, data=df_long_cbcl, 
                  boxwex=0.5, ylab="Canonical Loadings", outline=FALSE,
                  col=c("#A0D568","#F7EA48","#4FC1E8"), 
                  border=c("#1A9E37", "#F1B51A", "#3659A0"),
                  xaxt="n", alpha=0.7,ylim = c(-1, 1))
# Add the grey vertical lines
for(i in seq(0.5 , 24, 3)){ 
  abline(v=i,lty=1, col="grey")
}
abline(h=0, col="grey", lty=2)
# To add the label of x axis
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 3)]
axis(1, font=2,
     at = seq(2, 24, 3), 
     labels = my_names , 
     tick=FALSE , cex=0.3)
# Add a legend
legend(23.5,-0.2, legend = c("CV1", "CV2", "CV3"), 
       col=c("#A0D568","#F7EA48","#4FC1E8"),
       border=c("#1A9E37", "#F1B51A", "#3659A0"),
       pch = 15,  pt.cex = 2, cex = 0.9,  horiz = F, inset = c(0.1, 0.1))
dev.off()



######################################################
####### visualization of the CIs of brain loadings ####
######################################################
# training set
n_boot <- 1000

boot_brain <- lapply(1:n_boot, function(i) {
  boot_abcd[[i]]$brain_loading[,cv_reorder_boot_abcd[[i]]$idx_reorder]})

brain_sign <- sign(res.abcd.7$u)
max1 <- which(rank(-abs(res.abcd.7$u[,1])) == 1)
sign1 <- sign(res.abcd.7$u[max1,1])
max2 <- which(rank(-abs(res.abcd.7$u[,2])) == 1)
sign2 <- sign(res.abcd.7$u[max2,2])
max3 <- which(rank(-abs(res.abcd.7$u[,3])) == 1)
sign3 <- sign(res.abcd.7$u[max3,3])

sign <- c(sign1,sign2,sign3)

df_brain <- lapply(1:3, function(i) {
  brain_cv <- do.call(cbind, lapply(boot_brain, function(x){x[, i]}))
  brain_cv <- t((brain_cv))
  brain_cv <- as.data.frame(brain_cv)
  # adjust the signs 
  cvmean <- colMeans(abs(brain_cv))
  idx <- which(rank(-cvmean) == 1) # find the most important one
  
  for (j in 1:nrow(brain_cv)){
    if(sign(brain_cv[j, idx]) != sign[i]){
      brain_cv[j,] <- -brain_cv[j, ]
    }
  }
  names(brain_cv) <- paste0("PC",1:100)
  brain_cv$CV <- paste0("CV", i)
  return(brain_cv)})


df_brain <- do.call(rbind,df_brain)
names(df_brain)

# select the important PCs2
cvmean1 <- colMeans(abs(df_brain[df_brain$CV == "CV1", -ncol(df_brain)]))
cvmean2 <- colMeans(abs(df_brain[df_brain$CV == "CV2", -ncol(df_brain)]))
cvmean3 <- colMeans(abs(df_brain[df_brain$CV == "CV3", -ncol(df_brain)]))
cvmean1[order(-cvmean1)][1:10]
cvmean2[order(-cvmean2)][1:10]
cvmean3[order(-cvmean3)][1:10]
idx1 <- which(rank(-cvmean1) %in% 1:10)
idx2 <- which(rank(-cvmean2) %in% 1:10)
idx3 <- which(rank(-cvmean3) %in% 1:10)

# manually selected
PC_select <- c(1,2,3,4,10,17,26,44,86)
## final: 1,2,3,4,10,15,20,29

df_brain1 <- df_brain[, c(PC_select,ncol(df_brain))]
PCs <- paste0("PC", PC_select)

df_long_brain <- gather(df_brain1, PCs, loadings, PC1:PC86, factor_key = TRUE)


############ combine three CVs together
# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND treatment :
pdf(file="brain_loadings_boot_newsign.pdf", width=13, height=5)
myplot <- boxplot(loadings ~ CV*PCs, data=df_long_brain, 
                  boxwex=0.5, ylab="Canonical Loadings", outline=FALSE,
                  col=c("#A0D568","#F7EA48","#4FC1E8"),
                  border=c("#1A9E37", "#F1B51A", "#3659A0"),
                  xaxt="n", alpha=0.7,ylim = c(-1, 1))
# Add the grey vertical lines
for(i in seq(0.5, 27, 3)){ 
  abline(v=i,lty=1, col="grey")
}
abline(h=0, col="grey", lty=2)
# To add the label of x axis
my_names <- sapply(strsplit(myplot$names, '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 3)]
axis(1, font=2,
     at = seq(2, 27, 3), 
     labels = my_names , 
     tick=FALSE , cex=0.3)
# Add a legend
legend(26.5,0.9, legend = c("CV1", "CV2", "CV3"), 
       col=c("#A0D568","#F7EA48","#4FC1E8"),
       pch = 15,  pt.cex = 2, cex = 0.9,  horiz = F, inset = c(0.1, 0.1))
dev.off()



########################################################
############## Covariance explained  ###################
########################################################
vardf <- VarianceExplain(brain_train, cbcl_train, res.abcd.7, 8) 
colnames(vardf) <- c("Canonical_Variates", "Covariance_Explained")
vardf$`Canonical_Variates` <- paste0("CV",1:8)
vardf$Canonical_Variates <- as.factor(vardf$Canonical_Variates)
ggplot(vardf, aes(x=Canonical_Variates, y=Covariance_Explained, group=Canonical_Variates))+
  geom_point(aes(shape=Canonical_Variates,color=Canonical_Variates), size = 7, alpha = 0.9) + 
  theme_bw() + theme(legend.position="none",
                     axis.text.x = element_text(size = 12)) + 
  scale_shape_manual(values=c(17,16,15,18,20,20,8,8)) +
  scale_color_manual(values=c("#A0D568","#F7EA48","#4FC1E8","#f08080","#cc99c9","#ffce54","grey","grey")) 

ggsave("covariance_explained_trainset.pdf",width=8,height=4)


### permutation test: an example
perm_abcd7 <- permutation_test(cbcl_train, brain_train, nperm=999, 0.5, 0.5, 8, res.abcd.7$cors)


########################################################
############## permutation test visualization  #########
########################################################

# direct mapping
cor.abcd.test7 <- test_project_weights(brain_test, cbcl_test, res.abcd.7, 8)
cor.abcd.test7
perm_abcdtest7 <- permutation_test_testset(cbcl_test, brain_test,nperm=1999, 
                                           res.abcd.7,cor.abcd.test7)
#"#A0D568","#F7EA48","#4FC1E8"

p.adjust(perm_abcdtest7$pval.perm)

df.cor.2 <- data.frame(cor = perm_abcdtest7$cor.perm[, 2])
df.cor.1 <- data.frame(cor = perm_abcdtest7$cor.perm[, 1])
df.cor.3 <- data.frame(cor = perm_abcdtest7$cor.perm[, 3])

pdf("permutation_test_cv3.pdf", width=5, height=3)
ggplot(df.cor.3, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "#4FC1E8", color = "black", alpha = 0.8)+
  geom_vline(xintercept=0.075, linetype="dashed", color = "red", size = 1) +
  ggtitle("Permutation correlations of the 3rd canonical variate") + 
  labs(x = "canonical correlations") + 
  #annotate(geom="text", x=0.0, y=200, label="r = 0.13", size = , color="black") + 
  xlim(c(-0.2, 0.2)) + annotate(geom="text", x=0.15, y=180, label="P(FDR) < 0.01", 
                                fontface = 'italic', size = 5, color="black") + 
  theme_bw()
dev.off()



###############################################################
############## Visualization of correlations across 10 splits 
###############################################################
####### reorder the canonical variates

rs_train_test_abcd <- readRDS("rs_train_test_abcd.rds")
# the reference split
load_std <- res.abcd.7$v

new_cv_reorder  <- lapply(seq_along(rs_train_test_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:6) {
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:6, function(z) {cor(load_std[, i], rs_train_test_abcd[[j]]$abcd.train$v[, z])})
    if(max(abs(cor_cv)) < 0.5){ # sometimes the correlation is too low
      idx  <- i
    } else {
      idx <- which(abs(cor_cv) == max(abs(cor_cv)))
    }
    
    if(!idx %in% idx_v) {
      idx_v <- c(idx_v, idx) # sometimes it will overlap
    } else {
      idx_v <- c(idx_v, i) 
    }
  }
  idx_v[duplicated(idx_v)] <- c(1:6)[!c(1:6) %in% idx_v]
  
  loadMat <- rs_train_test_abcd[[j]]$abcd.train$v[, idx_v]
  return(list(loadMat=loadMat, idx_reorder=idx_v))
  
})

c <- lapply(new_cv_reorder, function(x) {x$loadMat})
c_mean <- do.call(cbind, lapply(1:6, 
                                function(i) {
                                  rowMeans(abs(do.call(cbind, 
                                                       lapply(c, function(x) {x[, i]})
                                  )))}))

genr_cor <- do.call(cbind,lapply(1:10, function(x) {
  i <- new_cv_reorder[[x]]$idx_reorder
  rs_train_test_abcd[[x]]$res.genr[i]}))
rowMeans(abs(genr_cor))


genr_cor_perm <- do.call(cbind,lapply(1:10, function(x) {
  i <- new_cv_reorder[[x]]$idx_reorder
  perm <- rs_train_test_abcd[[x]]$genr.perm[i]
  p.adjust(perm, method = "fdr")
  }))


##### brain average

b <- lapply(seq_along(new_cv_reorder), function(x) {
  loadMat <- rs_train_test_abcd[[x]]$abcd.train$u[, new_cv_reorder[[x]]$idx_reorder]})

brain_mean <- do.call(cbind, lapply(1:6, 
                                function(i) {
                                  rowMeans(abs(do.call(cbind, lapply(b, function(x) {x[, i]})
                                  )))}))
##### calculate the CV scores
brain_cvscores_average10 <- scale(brain_whole) %*% brain_mean
dim(brain_cvscores_average10)
saveRDS(brain_cvscores_average10, "brain_cvscores_average10.rds")



# reorder the correlations
traincor <- lapply(rs_train_test_abcd, function(x) {x$abcd.train$cors})
traincor <- do.call(rbind,traincor)
traincor_reorder <- lapply(1:10, function(i) {traincor[i, new_cv_reorder[[i]]$idx_reorder]})
traincor_reorder <- do.call(rbind,traincor_reorder)
rownames(traincor_reorder) <- paste0("split",1:10)
colnames(traincor_reorder) <- paste0("CV",1:6)

meantrain <- colMeans(traincor_reorder)
sdtrain <- apply(traincor_reorder, 2, sd)


testcor <- lapply(rs_train_test_abcd, function(x) {x$abcd.test})
testcor <- do.call(rbind,testcor)
testcor_reorder <- lapply(1:10, function(i) {testcor[i, new_cv_reorder[[i]]$idx_reorder]})
testcor_reorder <- do.call(rbind,testcor_reorder)
rownames(testcor_reorder) <- paste0("split",1:10)
colnames(testcor_reorder) <- paste0("CV",1:6)
meantest <- colMeans(testcor_reorder)
sdtest <- apply(testcor_reorder, 2, sd)

fig_10splits_syn <- data.frame(Train_Test = rep(c("Training set","Test set"),each=3), meancor = c(meantrain[1:3], meantest[1:3]),
                               sdcor = c(sdtrain[1:3], sdtest[1:3]), CV=rep(paste0("CV",1:3),2))


p <- ggplot(fig_10splits_syn, aes(x=CV, y=meancor, color=Train_Test, group=Train_Test)) + 
  scale_color_manual(name="Training Test Set", values = c("#fbb172", "#E64659")) + 
  geom_errorbar(aes(ymin=meancor-sdcor, ymax=meancor+sdcor), width=.2, size=2) +
  geom_point(size=5) + ylim(c(0, 0.25)) + 
  theme_bw() + labs(y = "canonical correlations", x = "canonical variates", size=5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text = element_text(size=13))
p

ggsave("train_test_variance.pdf", width = 8.5, height = 6)


####### extract average brain loadings

loadMat <- lapply(1:10, function(i) rs_train_test_abcd[[i]]$abcd.train$v[, idx_v])




###############################################################
############## Visualization of CBCL syndrome loadings  #######
###############################################################

loadings_cv <- lapply(new_cv_reorder, function(x) {x$loadMat})

# calculate the mean of across the splits
temp <- sapply(1:5, function(i) {
  rowMeans(sapply(seq_along(loadings_cv), function(j) {abs(loadings_cv[[j]][, i])}))
})

df_temp <- as.data.frame(temp)


##################################
######### circular radar plot
##################################
temp <- c_mean
df_temp1 <- as.data.frame(temp)
df_temp1 <- t(df_temp1)
df_temp2 <- rbind(rep(1,8), df_temp1)
colnames(df_temp2) <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule breaking","Aggression")

df_temp3 <- as.data.frame(cbind(group=c("1","CV1", "CV2", "CV3", "CV4", "CV5", "CV6"),df_temp2))
df_temp3[,2:9] <- lapply(df_temp3[,2:9], function(x){as.numeric(as.character(x))})
df_temp3 <- df_temp3[1:4,]

#pdf("radarfigure_train_abcd.pdf", width=12, height=8)
colors_border_6 <- c("#fcfcfb","#A0D568","#F7EA48","#4FC1E8","#F08080","#CC99C9","#FFCE54")
colors_in_6 <- c("#fcfcfb","#A0D568","#F7EA48","#4FC1E8","#F08080","#CC99C9","#FFCE54")

p2 <- ggradar2(df_temp3, base.size=5, webtype = 'lux', grid.min = 0,
         grid.max = 1,label.gridline.mid = TRUE, group.colours = colors_border_6,
         group.fill.colours = colors_in_6,label.centre.y=FALSE,
         gridline.mid.colour="grey", grid.label.size = 0,
         gridline.max.linetype = "solid",polygonfill.transparency=0.5,
         background.circle.transparency=0.1,axis.label.size=5.5,
         group.line.width=1.5,group.point.size=3,plot.legend = TRUE) 
p2
ggsave("radarfigure_test_abcd.png", width = 12, height = 6)
dev.off()




############################################
######### calculate the average correlations 
#############################################
rs_train_test_abcd <- readRDS("rs_train_test_abcd.rds")
# training set
colMeans(do.call(rbind, lapply(1:10, function(x) rs_train_test_abcd[[x]]$abcd.train$cors)))
# test set
colMeans(abs(do.call(rbind, lapply(1:10, function(x) rs_train_test_abcd[[x]]$res.genr))))

lapply(1:10, function(x) rs_train_test_abcd[[1]]$abcd.test.perm)



############################################
######### clustering 
#############################################

scca_brain <- readRDS("brain_cvscores_average10.rds")
cca_brain_1 <- scca_brain[, 1:3]

# 1. hiearachical clustering
d <- dist(cca_brain_1, method ="euclidean")
res.hc <- hclust(d, method = "ward.D")
clusters <- cutree(res.hc, k = 2)

df.clu.figure <- data.frame(CV1 = unlist(cca_brain_1[, 1]), 
                            CV2 = unlist(cca_brain_1[, 2]),
                            CV3 = unlist(cca_brain_1[, 3]))
df.clu.figure$clusters <- clusters
ggplot(df.clu.figure) + geom_point(aes(x = CV1, y = CV2), colour= clusters) +
  theme_bw()

# Dendrogram
d <- as.dendrogram(res.hc)
d <- d %>% color_branches(k=2) %>% color_labels
plot(d, main = "Dendrogram for hierarchical clustering")


# decide how many clusters
cca_rs_data <- cca_brain_1
hcfit_ch <- NbClust(cca_rs_data, method="ward.D",
                    min.nc = 1, max.nc = 6, index = "ch")
hcfit_sl <- NbClust(cca_rs_data, method="ward.D",
                    min.nc = 1, max.nc = 6, index = "silhouette")


fit_ch <- data.frame(number = 1:6, values = hcfit_ch$All.index)
fit_sl <- data.frame(number = 1:6, values = hcfit_sl$All.index)

# figures
ch <- ggplot(fit_ch, aes(number, values))+ 
  geom_line(color = "steelblue2")+ geom_point(color="navyblue" )+
  labs( x= "Number of clusters", y="variance ratio criterion") +
  ggtitle("Variance ratio criterion (Calinski-Harabasz index)") +
  theme_bw()
sl <- ggplot(fit_sl, aes(number, values))+ 
  geom_line(color = "hotpink1")+ geom_point(color="hotpink3" )+
  labs( x= "Number of clusters", y="Silhouette") +
  ggtitle("Silhoutte") +
  theme_bw()

grid.arrange(ch, sl, ncol=2)

# 2. Fuzzy C means clustre analysislibrary(cluster)

cluster_fuzzy_train <- cmeans(cca_brain_1 , centers=2, iter.max = 100, verbose = FALSE,
                              dist = "euclidean", method = "cmeans", m = 2,
                              rate.par = NULL, weights = 1, control = list())
cluster.df <- data.frame(Hclu=clusters, Fclu_1=cluster_fuzzy_train$membership[,1], 
                         Fclu_2 = cluster_fuzzy_train$membership[,2])
dim(cluster.df)


clu_stable <- cluster.df %>% filter(Fclu_1 > 0.7 | Fclu_2 > 0.7)
dim(clu_stable)

idx_clu <- which(cluster.df$Fclu_1 > 0.7 | cluster.df$Fclu_2 > 0.7)

RS_stable <- data.frame(CV1 = unlist(cca_brain_1[-idx,1]), CV2 = unlist(cca_brain_1[-idx,2]),
                        CV3 = unlist(cca_brain_1[-idx,3]))
RS_stable$clusters <- as.factor(clusters[-idx])

idx <- which(cca_brain_1[,2]== min(cca_brain_1[,2]))




plot_ly(RS_stable, x = ~CV1, y = ~CV2, z = ~CV3, 
        color = ~clusters, marker = list(size=3.5), colors = c("skyblue2", "gold1"))

dim(cca_rs_data)
# significance test
sigma <- cov(cca_rs_data)
mu <- colMeans(cca_rs_data)
real_CI <- cluster_test(cca_rs_data)

cl <- makePSOCKcluster(8)
registerDoParallel(cl)
n_sims <- 99
null <- foreach(i = 1:n_sims, .packages=c("NbClust","MASS")) %dopar% {
  rand_sample <- mvrnorm(n=nrow(cca_rs_data), mu=mu, Sigma=sigma)
  res <- cluster_test(rand_sample)
}

stopCluster(cl)



cluster_test <- function(cca_data){
  #ugly hack, because i don't know how to prevent this library creating many plots
  hcfit <- NbClust(cca_data, method="ward.D", index="ch", min.nc=1, max.nc = 6)
  CH_index <- max(hcfit$All.index, na.rm = TRUE)
  hcfit <- NbClust(cca_data, method="ward.D", index="silhouette", min.nc=1, max.nc = 6)
  sil_index <- max(hcfit$All.index, na.rm = TRUE)
  return(c("CH"=CH_index, "Silhouette"=sil_index))
}






null_CI <- as.data.frame(do.call(rbind, null))

rank_cv1 <- sum(real_CI[1] < null_CI[,1]) + 1
pval_cv1 <- rank_cv1 / (n_sims+1)
rank_cv2 <- sum(real_CI[2] < null_CI[,2]) + 1
pval_cv2 <- rank_cv2 / (n_sims+1)
t(t((c("p.val variance ratio"=pval_cv1, "p.val Silhouette"=pval_cv2))))


df.ch <- data.frame(cor = null_CI[,1])
df.sl <- data.frame(cor = null_CI[,2])

ch1 <- ggplot(df.ch, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=20, fill = "gold1", color = "black", alpha = 0.5)+ 
  geom_vline(xintercept=real_CI[1], linetype="dashed", color = "red", size = 1) +
  ggtitle("Variance ratio criterion null") + 
  labs(x = "null distribution") + 
  annotate(geom="text", x=1300, y=60, label="p < 0.001", fontface = 'italic', size = 5, color="black") + 
  theme_bw()
sl2 <- ggplot(df.sl, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "lightblue", color = "black", alpha = 0.5)+
  geom_vline(xintercept=real_CI[2], linetype="dashed", color = "red", size = 1) + 
  annotate(geom="text", x=0.6, y=145, label="p < 0.001", fontface = 'italic', size = 5, color="black") +
  ggtitle("Silhouette null") + 
  labs(x = "null distribution") + 
  theme_bw()

grid.arrange(ch1, sl2, ncol =2)




###############################################
######### clustering: demogrpahic distributions 
###############################################
setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC")

all_demo <- readRDS("all_final_noNA_incidental_6529.rds")
all_demo$clusters <- clusters

demo_cbcl <- all_demo[,c("sex","race_ethnicity","parental_education",
                        "cbcl_scr_syn_internal_r","cbcl_scr_syn_external_r",
                        "cbcl_scr_syn_totprob_r","clusters")]
summary(demo_cbcl[demo_cbcl$clusters == 1,])
lapply(demo_cbcl[demo_cbcl$clusters == 1,4:6], sd)





