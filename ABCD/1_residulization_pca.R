###################### feature selection ###############
library(doParallel)
library(permute)


############### functions
make_brain_features <- function(subid){
    
    conMatDir <- ('/gpfs/work2/0/einf1049/scratch/bxu/ABCD_data_filtered/')
    cl <- makePSOCKcluster(32)
    registerDoParallel(cl)
    
    feature_list <- foreach::foreach(i = 1:length(subid), .packages=c("doParallel", "foreach")) %dopar% {
        subID <- subid[i]
        conMatCsv <- file.path(conMatDir, paste(subID, 'corMat_Gordon_fsSC_filtered', sep='_'))
        conMatCsv <- paste(conMatCsv, 'csv', sep='.')
        conMat <- as.matrix(read.table(conMatCsv))
        feature_column <- conMat[upper.tri(conMat, diag = FALSE)]
    }
    
    stopCluster(cl)
    
    feature_abcd <- do.call(rbind,feature_list)
    feature_abcd <- as.data.frame(feature_abcd)
    colnames(feature_abcd) <- paste0("F_", 1:ncol(feature_abcd))
    
    return(feature_abcd)
}


weighted_pca <- function(cbcl,brain,n){
    rank <- rank(-rowSums(cbcl))
    weights <- log(nrow(cbcl)) - log(rank)
    weights <- weights / sum(weights)
    # center the original data
    feature_brain_centered <- scale(brain, scale = FALSE)
    # scale the columns of the feature matrix
    feature_centered_weighted <- feature_brain_centered * replicate(ncol(feature_brain_centered), weights)
    # call PCA on the weighted feature matrix
    pca.weighted <- prcomp(feature_centered_weighted)
    # the rotation transformation: take the first 100 eigenvectors
    rotation <- pca.weighted$rotation[,1:n]
    # the data with the reduced dimensionality
    feature_brain_reduced <- feature_brain_centered %*% rotation 
    
    return(list(brain_train_reduced=feature_brain_reduced, rotation = pca.weighted$rotation))
}


residualization <- function(brain,confounders){
    cl <- makePSOCKcluster(32)
    registerDoParallel(cl)
    
    residual_list <- foreach::foreach(i = seq_along(brain), .packages=c("doParallel", "foreach")) %dopar% {
        out <- residuals(lm(brain[, i] ~ confounders$age + confounders$sex + 
                          confounders$race_ethnicity + confounders$site + confounders$parental_education, na.action=na.exclude))
    }
    
    stopCluster(cl)
    
    brain_residual <- do.call(cbind, residual_list)
    return(brain_residual)
}



############### read the data

train <- readRDS("all_final_train.rds")
test <- readRDS("all_final_test.rds") 

train_test_split <- lapply(1:10, function(i) {

    train0 <- train[[i]]
    subid_train <- train0$participant_id
    subid_train <- as.character(subid_train)
    
    test0 <- test[[i]]
    subid_test <- test0$participant_id
    subid_test <- as.character(subid_test)
    
    ########### data extraction
    brain_train <- make_brain_features(subid_train)
    brain_test <- make_brain_features(subid_test)
    
    cbcl_train <- train0[, c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                           "cbcl_scr_syn_social_r", "cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                           "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")]
    cbcl_test <- test0[,c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                           "cbcl_scr_syn_social_r", "cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                           "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")]
                           
    confounders_train <- train0[,c("age","sex","race_ethnicity","site","parental_education")]
    confounders_test <- test0[,c("age","sex","race_ethnicity","site","parental_education")]
        
    ############ residualization
    brain_train_residual <- residualization(brain_train, confounders_train)
    brain_test_residual <- residualization(brain_test, confounders_test)
    
    ############ weighted PCA
    pca.weighted.rotation <- weighted_pca(cbcl_train, brain_train_residual, 100)
    
    out <- list(pca_train=pca.weighted.rotation, cbcl_train=cbcl_train, cbcl_test=cbcl_test, brain_test=brain_test_residual)
})


saveRDS(train_test_split, "train_test_split.rds")


