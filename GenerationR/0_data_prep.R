#####################################################################
################## Generation R: data preparation ###################
#####################################################################

library(stringr)
library(data.table)
library(base)
library(dplyr)

##########################################
###### rs-fMRI data QC ###################
##########################################


######## load all the inclusion data
setwd("V:/medewerkers/051950 Xu, B/PhD projects/Biotypes of ADHD/New inclusion criteria")
load("dataf9_aug20.RData")
load("motionf9_aug20.RData")
regis <- fread("registration_9_13.csv")
 
motion <- motionf9 %>% filter(mean_rms9 <= 0.25 & rms_vols_bin9 == 0 & !is.na(mean_rms9))
registration <- regis %>% filter((is.na(rsfmri_ants_reg_F09_1) | rsfmri_ants_reg_F09_1 == 8 |
                                   rsfmri_ants_reg_F09_1 == 2 | rsfmri_ants_reg_F09_1 == 7) &
                                   (is.na(rsfmri_ants_reg_F09_2) | rsfmri_ants_reg_F09_2 == 8 |
                                   rsfmri_ants_reg_F09_2 == 2 | rsfmri_ants_reg_F09_2 == 7 ))
incomplete_vol <- dataf9$idc[dataf9$dataf9 == 1]
 
rs_qc_pass_id <- incomplete_vol[incomplete_vol %in% intersect(motion$idc, registration$IDC)]
length(rs_qc_pass_id)
 
##########################################
###### merge CBCL ########################
##########################################
cbcl_total <- fread("CHILDCBCL9_incl_Tscores_20201111.csv")
cbcl_syndrom <- cbcl_total[, c("IDC","sum_anx_9m","sum_wit_9m","sum_som_9m","sum_sop_9m",
                               "sum_tho_9m","sum_att_9m","sum_rul_9m","sum_agg_9m")]
 
qc_cbcl <- cbcl_syndrom %>% filter(IDC %in% rs_qc_pass_id)
qc_cbcl[qc_cbcl == 888] <- NA
qc_cbcl[qc_cbcl == 999] <- NA
 
# check the NAs
table(colSums(is.na(qc_cbcl)))
 
# index of the participants with missingness in CBCL > 25%, 8 * 0.25 = 2
idx_cbcl_NA <- which(rowSums(is.na(qc_cbcl)) > 2 )
 
qc_cbcl_final <- qc_cbcl[-idx_cbcl_NA, ]
dim(qc_cbcl_final)
 
 
##########################################
###### remove one of twins or siblings ###
##########################################
 
demo <- fread("CHILD-ALLGENERALDATA_12112020.csv")
demo_qc_cbcl <- merge(qc_cbcl_final,demo[,c("IDC","MOTHER")],by="IDC")
dim(demo_qc_cbcl)
demo_qc_cbcl_final <- demo_qc_cbcl[!duplicated(demo_qc_cbcl$MOTHER), ]
dim(demo_qc_cbcl_final)

##########################################
###### add covariates ####################
##########################################

MRI_core <- readRDS("genr_mri_core_data_20220311.rds") # mri related demographic data
names(MRI_core)[1] <- "IDC"
databrain <- fread("F9_MRI_Freesurfer_global_2016_12_09.csv") # freesurfer data
dataE <- fread("SESINTAKE_17072015.csv") # intake education

####### merge the specific data need for imputation
### mri demographic info
df_1 <- merge(demo_qc_cbcl_final, MRI_core[, c("IDC","age_child_mri_f09")], by="IDC")
dim(df_1)
### other SES info
df_2 <- merge(df_1, demo[, c("IDC", "ETHNINFv3", "GENDER", "AGE_M_v2", "INCOME5", "GESTBIR")], by="IDC")
### total brain volume
df_3 <- merge(df_2, databrain[, c("IDC", "knicr_tbv_f9")], by="IDC")


############## special cases one: maternal education
##### This dataset does not have the IDs for children, just for mothers
##### and maternal education was assessed in multiple times, so I will merge 
##### all the maternal education variables from 5 years old, 36 months and prenatal

######### maternal education
maternal_edu <- merge(dataAll[, c("IDC","IDM","EDUCM5","EDUCM3")], dataE[,c("IDM", "EDUCM")], by="IDM")
maternal_edu[maternal_edu == 888] <- NA
maternal_edu[maternal_edu == 999] <- NA
# when all three time points are NA, assign NA, otherwise choose the highest one
edu_m <- apply(maternal_edu[,3:5], 1, function(x) {ifelse(sum(is.na(x))==3, NA, max(x[!is.na(x)]))})
edum_final <- cbind(maternal_edu, edu_m)

df_4 <- merge(df_3, edum_final[, c("IDC", "edu_m")], by="IDC")
dim(df_4)

############## special case two: re-level ethnicity
# Dutch
# non-dutch european: turkish, Surinamese_Hindustani, european
# Other: Moroccan, Indonesian, Cape Verdian, Dutch Antilles, Surinamese,
# Surinamese-Creole, Surinamese_unspecified,African, American western, American non Western, Asian western
# Asian non western, Oceanie. 

df_4[df_4 == 888] <- NA
df_4[df_4 == 999] <- NA

df_4$ETHNINFv3 <- factor(df_4$ETHNINFv3)
levels(df_4$ETHNINFv3) <- list( "dutch" = 1,
                                "european" = c(7,9,700),
                                "other" = c(2,3,4,5,6,8,10,200,300,400,500,600,800))
table(df_4$ETHNINFv3)
df_final <- df_4
dim(df_final)

##################################################
################### EM imputation ################
##################################################
library(Amelia)
library(dplyr)

###### merge data and adjust the data struture

# read the brain PCs
PCs_genr <- readRDS("PCs_genr.rds")
dim(PCs_genr)
# use the first 10 brain PCs used for imputation
df_final <- cbind(df_final, PCs_genr[,1:10])

df_final_imp <- df_final %>% mutate(age = scale(age_child_mri_f09), age_pre=scale(AGE_M_v2),
                                    gender=GENDER, ethnicity=factor(ETHNINFv3), income=INCOME5,
                                    gestbir = scale(GESTBIR), maternal_edu=edu_m,
                                    totalbrain=scale(knicr_tbv_f9))
names(df_final_imp)
# remove the original duplicated variables
df_final_imp <- df_final_imp[, !names(df_final_imp) %in% c("age_child_mri_f09","AGE_M_v2","GENDER","ETHNINFv3",
                                                           "INCOME5","GESTBIR","edu_m","knicr_tbv_f9","IDC")]
###### single imputation
all_toimp <- amelia(df_final_imp, noms=c("ethnicity"), ords=c("maternal_edu", "income"), m=1, boot.type = "none")

# check the output
lapply(all_toimp$imputations$imp1, table)

final_imp <- cbind(IDC=df_final$IDC, all_toimp$imputations$imp1)

final_imp$maternal_edu <- factor(final_imp$maternal_edu)
# the levels should be changed after the imputation
levels(final_imp$maternal_edu) <- list( "low" = 0:1,
                                        "medium" = 2:3,
                                        "high" = 4:5)

saveRDS(final_imp, "final_imp_cbcl_covariates.rds")
