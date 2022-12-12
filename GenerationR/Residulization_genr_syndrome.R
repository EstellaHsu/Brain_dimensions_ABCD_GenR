###################### feature selection ###############
setwd("/home/r051950/Gordon_newdata/Gordon_AROMAOnly/")
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
