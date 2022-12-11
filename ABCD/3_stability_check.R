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


#############################################################################
###### bootstrap of one split: stability check of cbcl and the brain loaidngs
#############################################################################

######## read the data: for example split 7
rs_example <- readRDS("train_test_split.rds")
rs_example <- rs_example[[7]]
brain_train <- rs_example$pca_train$brain_train_reduced
cbcl_train <- rs_example$cbcl_train
brain_test <- as.matrix(rs_example$brain_test) %*% rs_example$pca_train$rotation[, 1:100]
cbcl_test <- rs_example$cbcl_test 
# penalty parameters are based on grid search in the second SCCA step. 
res.abcd.example <- CCA(x=brain_train, z=cbcl_train, penaltyx = 0.5, penaltyz = 0.5, typex="standard", typez="standard",
                  niter = 20, K=8)
res.abcd.example


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
load_std <- res.abcd.example$v
cv_reorder_boot_abcd  <- lapply(seq_along(boot_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:5) { # we only care about the first 5 CVs, as the covariance explained is very small for other CVs
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:5, function(z) {cor(load_std[, i], boot_abcd[[j]]$cbcl_loading[, z])})
    if(max(abs(cor_cv)) < 0.5){ # sometimes the correlation is too low
      idx  <- i
    } else {
      idx <- which(abs(cor_cv) == max(abs(cor_cv)))
    }
    
    if(!idx %in% idx_v) {# sometimes it will overlap
      idx_v <- c(idx_v, idx) 
    } else {
      idx_v <- c(idx_v, i) 
    }
  }
  
  idx_v[duplicated(idx_v)] <- c(1:5)[!c(1:5) %in% idx_v]
  list(idx_reorder=idx_v)  
})


##### first three canonical correlations in training sets
boot_cor_train <- lapply(boot_abcd, function(x) {x$cors_train[1:3]})
boot_cor_train <- do.call(rbind,boot_cor_train)
boot_cor_train_reorder <- lapply(1:n_boot, function(i) {boot_cor_train[i, cv_reorder_boot_abcd[[i]]$idx_reorder]})
boot_cor_train_reorder <- do.call(rbind, boot_cor_train_reorder)
colnames(boot_cor_train_reorder) <- paste0("CV",1:3)
meantrain_boot <- colMeans(boot_cor_train_reorder)
sdtrain_boot <- apply(boot_cor_train_reorder,2,sd)

##### first three canonical correlations in test sets
boot_cor_test <- lapply(boot_abcd, function(x) {x$cors_test[1:3]})
boot_cor_test <- do.call(rbind,boot_cor_test)
boot_cor_test_reorder <- lapply(1:n_boot, function(i) {boot_cor_test[i, cv_reorder_boot_abcd[[i]]$idx_reorder]})
boot_cor_test_reorder <- do.call(rbind,boot_cor_test_reorder)
colnames(boot_cor_test_reorder) <- paste0("CV",1:3)
meantest_boot <- colMeans(boot_cor_test_reorder)
sdtest_boot <- apply(boot_cor_test, 2, sd)


#################################################
########### violin plot with error bars
#################################################
#train set
df1 <- as.data.frame(boot_cor_train_reorder[,1:3])
names(df1) <- paste0("CV", 1:3)
df1$traintest <- "Train"
#test set
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
ggsave("train_test_boot_violin.pdf", width = 8, height = 6)


######################################################
####### visualization of the CIs of cbcl loadings ####
######################################################

cbclnames <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule Breaking","Aggression")
# training set
n_boot <- 1000
boot_cbcl <- lapply(1:n_boot, function(i) {
  boot_abcd[[i]]$cbcl_loading[,cv_reorder_boot_abcd[[i]]$idx_reorder]})

cbcl_sign <- sign(res.abcd.example$v)

# cv1: attention, cv2: aggression & rule breaking, cv3:withdrawn

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


############ combine three CVs together and make the figures

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


#########################################################
####### visualization of the CIs of brain loadings ######
#########################################################
                   
# training set
n_boot <- 1000
                   
# extract brain loadings
boot_brain <- lapply(1:n_boot, function(i) {
  boot_abcd[[i]]$brain_loading[,cv_reorder_boot_abcd[[i]]$idx_reorder]})

                   
brain_sign <- sign(res.abcd.example$u)
# because there are 100 brain features, we just make sure the most important one get the correct sign
max1 <- which(rank(-abs(res.abcd.example$u[,1])) == 1)
sign1 <- sign(res.abcd.example$u[max1,1])
max2 <- which(rank(-abs(res.abcd.example$u[,2])) == 1)
sign2 <- sign(res.abcd.example$u[max2,2])
max3 <- which(rank(-abs(res.abcd.example$u[,3])) == 1)
sign3 <- sign(res.abcd.example$u[max3,3])

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

# select the important brain PCs
cvmean1 <- colMeans(abs(df_brain[df_brain$CV == "CV1", -ncol(df_brain)]))
cvmean2 <- colMeans(abs(df_brain[df_brain$CV == "CV2", -ncol(df_brain)]))
cvmean3 <- colMeans(abs(df_brain[df_brain$CV == "CV3", -ncol(df_brain)]))
cvmean1[order(-cvmean1)][1:10]
cvmean2[order(-cvmean2)][1:10]
cvmean3[order(-cvmean3)][1:10]

# manually select the top 3 PCs for each CV
PC_select <- c(1,2,3,4,10,17,26,44,86)

df_brain1 <- df_brain[, c(PC_select,ncol(df_brain))]
PCs <- paste0("PC", PC_select)

df_brain1 <- df_brain[, c(PC_select,ncol(df_brain))]

df_long_brain <- gather(df_brain1, PC_select, loadings, PC1:PC86, factor_key = TRUE)

############ combine three CVs together

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





