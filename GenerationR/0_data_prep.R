#####################################################################
################## Generation R: data preparation ###################
#####################################################################

library(stringr)
library(data.table)
library(base)
library(dplyr)

#################### data extraction ##############

MRI_core <- readRDS("genr_mri_core_data_20220311.rds") # mri related demographic data
names(MRI_core)[1] <- "IDC"
dataAll <- fread("CHILD-ALLGENERALDATA_12112020.csv") # other demographic data 
CBCL <- fread("CHILDCBCL9_incl_Tscores_20201111.csv") # CBCL scores
databrain <- fread("F9_MRI_Freesurfer_global_2016_12_09.csv") # freesurfer data
dataE <- fread("SESINTAKE_17072015.csv") # intake education

####### extract the specific data need for imputation

### demographic info
df_1 <- merge(cbcl_final, MRI_core[, c("IDC","age_child_mri_f09")], by="IDC")
dim(df_1)
### other SES info
df_2 <- merge(df_1, dataAll[, c("IDC", "ETHNINFv3", "GENDER", "AGE_M_v2", "INCOME5", "GESTBIR")], by="IDC")
### CBCL syndrome scale scores
df_3 <- merge(df_2, CBCL[, c("IDC", "sum_anx_9m","sum_wit_9m","sum_som_9m","sum_sop_9m",
                                             "sum_tho_9m","sum_att_9m","sum_rul_9m","sum_agg_9m")], by="IDC")
dim(df_3)
### total brain volume
df_4 <- merge(df_3, databrain[, c("IDC", "knicr_tbv_f9")], by="IDC")


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

df_5 <- merge(df_4, edum_final[, c("IDC", "edu_m")], by="IDC")
dim(df_5)

############## special case two: re-level ethnicity
# Dutch
# non-dutch european: turkish, Surinamese_Hindustani, european
# Other: Moroccan, Indonesian, Cape Verdian, Dutch Antilles, Surinamese,
# Surinamese-Creole, Surinamese_unspecified,African, American,western, American,non Western, Asian, western
# Asian, non western, Oceanie. 

df_5[df_5 == 888] <- NA
df_5[df_5 == 999] <- NA

df_5$ETHNINFv3 <- factor(df_5$ETHNINFv3)
levels(df_5$ETHNINFv3) <- list( "dutch" = 1,
                                "european" = c(7,9,700),
                                "other" = c(2,3,4,5,6,8,10,200,300,400,500,600,800))
table(df_5$ETHNINFv3)

df_final <- df_5
dim(df_final)

##################################################
################### EM imputation ################
##################################################
library(Amelia)
library(dplyr)


##############################################
###### merge data and adjust the data struture
setwd("V:/medewerkers/051950 Xu, B/PhD projects/TotalCBCL/GenR_Gordon_AROMAOnly")
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
# single imputation
all_toimp <- amelia(df_final_imp, noms=c("ethnicity"), ords=c("maternal_edu", "income"), m=1, boot.type = "none")

# check the output
lapply(all_toimp$imputations$imp1, table)

final_imp <- cbind(IDC=df_final$IDC, all_toimp$imputations$imp1)
final_imp$maternal_edu <- factor(final_imp$maternal_edu)
# the levels should be changed after the imputation
levels(final_imp$maternal_edu) <- list( "low" = 0:1,
                                        "medium" = 2:3,
                                        "high" = 4:5)


saveRDS(final_imp, "confounders_imp_cbclsyndrome_meduRelevel_aromaonly_notwinsib.rds")

