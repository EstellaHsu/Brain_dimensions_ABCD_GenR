###################### feature selection ###############
setwd("/gpfs/work2/0/einf1049/scratch/bxu/ABCD_replication_codes/finalReleaseQC_0817/")
library(doParallel)
library(permute)

brain <- readRDS("feature_abcd_finalReleaseQC_noNA.rds")

confounders <- readRDS("all_demo_qc_gene_cbclsyndromeRawscore_newQC.rds")

id <- intersect(brain$participant_id, confounders$participant_id)

brain_filter <-  brain[brain$participant_id %in% id, ]
brain_filter <- brain_filter[, -ncol(brain_filter)]

confounders <- confounders[confounders$participant_id %in% id, ]


# regress the feature matrix on age, gender, ethnicity, sites and SES (maternal education)
cl <- makePSOCKcluster(32)
registerDoParallel(cl)

residual_list <- foreach::foreach(i = seq_along(brain_filter), .packages=c("doParallel", "foreach")) %dopar% {
    out <- residuals(lm(brain_filter[, i] ~ confounders$age + confounders$sex + 
                      confounders$race_ethnicity + confounders$site + confounders$parental_education, na.action=na.exclude))
}

stopCluster(cl)

brain_residual <- do.call(cbind, residual_list)
brain_residual <- cbind(participant_id=id, brain_residual)

saveRDS(brain_residual, "brain_residual_abcd_cbclsyndrome_finalReleaseQC.rds")
saveRDS(id, "id_intersectBrainCbclsyndrome.rds")