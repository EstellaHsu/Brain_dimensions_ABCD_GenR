setwd("/gpfs/work2/0/einf1049/scratch/bxu/ABCD_replication_codes/finalReleaseQC_0817/")
library(doParallel)
library(permute)


behavior <- readRDS("all_demo_qc_gene_cbclsyndromeRawscore_newQC.rds")
subid <- behavior$participant_id
subid <- as.character(subid)

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
feature_abcd$participant_id <- subid
saveRDS(feature_abcd, "feature_abcd_finalReleaseQC_6816.rds")






