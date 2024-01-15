###############################################################################
#################### GenR: residualization & weighted PCA #####################
###############################################################################


library(doParallel)
library(permute)

#####################################
######## 1.read the brain data ######
#####################################
cbcl_demo <- readRDS("final_imp_cbcl_covariates.rds")
subid <- cbcl_demo$IDC
conMatDir <- ('PATH TO YOUR DATA')
 
cl <- makePSOCKcluster(12)
registerDoParallel(cl)
 
feature_list <- foreach::foreach(i = 1:length(subid), .packages=c("doParallel", "foreach")) %dopar% {
    subID <- paste('sub', subid[i], sep='-')
    conMatCsv <- file.path(conMatDir, paste0(subID, '_corMatZ.Gordon_fsSC.8pAroma.txt'))
    conMat <- as.matrix(read.table(conMatCsv))
    feature_column <- conMat[upper.tri(conMat, diag = FALSE)]
}
 
stopCluster(cl)
 
feature_genr <- do.call(rbind,feature_list)
feature_genr <- as.data.frame(feature_genr)
colnames(feature_genr) <- paste0("F_", 1:ncol(feature_genr))


#####################################
######## 2.residualization ##########
#####################################

# regress the feature matrix on age, gender, ethnicity and SES (maternal education)
resi_brain <- lapply(1:(ncol(brain)), function(i) {
   out <- residuals(lm(feature_genr[, i] ~ cbcl_demo$age + cbcl_demo$gender + 
                      cbcl_demo$ethnicity + cbcl_demo$maternal_edu, na.action=na.exclude))
})

brain_residual <- do.call(cbind, resi_brain)
colnames(brain_residual) <- paste0("F_", 1:ncol(brain_residual))

#####################################
######## 2.residualization ##########
#####################################
cbcl_genr <- cbcl_demo[, c("sum_anx_9m","sum_wit_9m","sum_som_9m","sum_sop_9m",
                           "sum_tho_9m","sum_att_9m","sum_rul_9m","sum_agg_9m")]
pca.weighted.genr <- weighted_pca(cbcl_genr, brain_residual, 100)
# extract the first 100 brain PCs
brain_genr <- pca.weighted.genr$brain_train_reduced






