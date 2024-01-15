####################################################################################
###################### Read the data, residualization, weighted PCA ################
####################################################################################

library(doParallel)
library(permute)

############### read the data (30 train-test splits)

train <- readRDS("all_final_train.rds")
test <- readRDS("all_final_test.rds") 

############### residualization and weighted PCA separately for ABCD training and ABCD test 

train_test_split <- lapply(1:30, function(i) {

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




