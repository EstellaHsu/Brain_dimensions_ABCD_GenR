###############################################################################
#################### GenR: residualization & weighted PCA #####################
###############################################################################



#####################################
######## 1.residualization ##########
#####################################

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












brain <- readRDS("feature_genr_aromaonly_notwinsib.rds")
brain <- brain[, -ncol(brain)]
confounders_imp <- readRDS("confounders_imp_cbcl_meduRelevel_aromaonly_notwinsib.rds")

# regress the feature matrix on age, gender, ethnicity and SES (maternal education)
resi_brain <- lapply(1:(ncol(brain)), function(i) {
   out <- residuals(lm(brain[, i] ~ confounders_imp$age + confounders_imp$gender + 
                      factor(confounders_imp$ethnicity) + factor(confounders_imp$maternal_edu), na.action=na.exclude))
})


brain_residual <- do.call(cbind, resi_brain)
colnames(brain_residual) <- paste0("F_", 1:ncol(brain_residual))

saveRDS(brain_residual, "feature_residual_aromaonly_cbclitem.rds")
