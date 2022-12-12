library("stringr")
library("data.table")
library(base)
library(dplyr)
#################### confounders extraction ##############


setwd("V:/medewerkers/051950 Xu, B/PhD projects/TotalCBCL/GenR_Gordon_AROMAOnly")

cbcl_final <- readRDS("cbcl_final.rds")
id <- readRDS("id_finalfinal_intersectBrainCBCL_removeTwinsSiblings.rds")
identical(as.character(id), as.character(cbcl_final$IDC))



################# remove the IDs in the covraties #####################
# extract all the data 
MRI_core <- readRDS("genr_mri_core_data_20220311.rds")
names(MRI_core)[1] <- "IDC"
dataAll <- fread("CHILD-ALLGENERALDATA_12112020.csv") # all demographic data 
CBCL <- fread("CHILDCBCL9_incl_Tscores_20201111.csv")
dataDrink <- fread("20141117_GEDRAGSGROEP_MaternalDrinking_UPDATE.csv")
databrain <- fread("F9_MRI_Freesurfer_global_2016_12_09.csv")
dataE <- fread("SESINTAKE_17072015.csv") # intake education



####### extract the data need for imputation
### demographic info
confounder_1 <- merge(cbcl_final, MRI_core[, c("IDC","age_child_mri_f09")], by="IDC")
dim(confounder_1)
### other SES info
confounder_2 <- merge(confounder_1, dataAll[, c("IDC", "ETHNINFv3", "GENDER", "AGE_M_v2", "INCOME5", "MAR_DICH5", "GESTBIR")], by="IDC")
### attention subscale
confounder_3 <- merge(confounder_2, CBCL[, c("IDC", "sum_anx_9m","sum_wit_9m","sum_som_9m","sum_sop_9m",
                                             "sum_tho_9m","sum_att_9m","sum_rul_9m","sum_agg_9m")], by="IDC")
dim(confounder_3)
### total brain volume
confounder_4 <- merge(confounder_3, databrain[, c("IDC", "knicr_tbv_f9")], by="IDC")
### ADHD medicine
#confounder_5 <- merge(confounder_4, dataADHD[, c("IDC", "ADHDmed")], by="IDC")


############## special cases one: maternal drink and maternal education
##### These two datasets don't have the ID for children, just for mother
##### and maternal education was assessed in multiple times, so I will merge 
##### all the maternal education variables from 5 years old, 36 months and prenatal

######### maternal education

maternal_edu <- merge(dataAll[, c("IDC","IDM","EDUCM5","EDUCM3")], dataE[,c("IDM", "EDUCM")], by="IDM")
maternal_edu[maternal_edu == 888] <- NA
maternal_edu[maternal_edu == 999] <- NA

edu_m <- apply(maternal_edu[,3:5], 1, function(x) {ifelse(sum(is.na(x))==3, NA, max(x[!is.na(x)]))})
edum_final <- cbind(maternal_edu, edu_m)

confounder_6 <- merge(confounder_4, edum_final[, c("IDC", "edu_m")], by="IDC")
dim(confounder_6)
######### maternal drink

#maternal_drink <- merge(dataAll[, c("IDC","IDM")], dataDrink[, c("IDM", "mdrink_updated")], by="IDM")
#confounder_7 <- merge(confounder_6, maternal_drink[, c("IDC", "mdrink_updated")], by="IDC")
#dim(confounder_7)

confounder_7 <- confounder_6

############## special case two: re-level ethnicity
# 1: Dutch
# 2: non-dutch european: turkish, Surinamese_Hindustani, european
# 3: Other: Moroccan, Indonesian, Cape Verdian, Dutch Antilles, Surinamese,
# Surinamese-Creole, Surinamese_unspecified,African, American,western, American,non Western, Asian, western
# Asian, non western, Oceanie. 
confounder_7[confounder_7 == 888] <- NA
confounder_7[confounder_7 == 999] <- NA

table(confounder_7$ETHNINFv3)
confounder_7$ETHNINFv3 <- factor(confounder_7$ETHNINFv3)
# the levels should be changed after the imputation
levels(confounder_7$ETHNINFv3) <- list( "dutch" = 1,
                                        "european" = c(7,9,700),
                                        "other" = c(2,3,4,5,6,8,10,200,300,400,500,600,800))
table(confounder_7$ETHNINFv3)



confounder_final <- confounder_7
dim(confounder_final)

saveRDS(confounder_final, "confounder_final.rds")
dim(confounder_final)
confounder_final <- readRDS("confounder_final.rds")


################### EM imputation ###############
library(Amelia)
library(dplyr)
#########################################################################
###################### confounders imputation ############################

setwd("V:/medewerkers/051950 Xu, B/PhD projects/TotalCBCL/GenR_Gordon_AROMAOnly")

PCs_genr <- readRDS("PCs_genr.rds")
dim(PCs_genr)
colSums(is.na(confounder_final))

confounder_final <- cbind(confounder_final, PCs_genr[,1:10])

table(confounder_final$INCOME5)
confounders <- confounder_final %>% mutate(age = scale(age_child_mri_f09), age_pre=scale(AGE_M_v2),
                                           gender=GENDER, ethnicity=factor(ETHNINFv3), income=INCOME5,
                                           marital_stat=MAR_DICH5, gestbir = scale(GESTBIR),
                                           maternal_edu=edu_m, 
                                           totalbrain=scale(knicr_tbv_f9))
names(confounders)
# remove the original duplicated variables
confounders <- confounders[, -c(121:127, 136:137)]
# remove ID and item-level data
df2 <- confounders[,-c(1:120)]
names(df2)

lapply(df2, table)



all_toimp <- amelia(df2, noms=c("ethnicity", "marital_stat"), 
                    ords=c("maternal_edu", "income"),
                    m=1, boot.type = "none")


lapply(all_toimp$imputations$imp1[, 22:26], table)


confounders_imp <- cbind(IDC=confounders$IDC, all_toimp$imputations$imp1)
dim(confounders_imp)


table(confounders_imp$maternal_edu)
confounders_imp$maternal_edu <- factor(confounders_imp$maternal_edu)
# the levels should be changed after the imputation
levels(confounders_imp$maternal_edu) <- list( "low" = 0:1,
                                             "medium" = 2:3,
                                             "high" = 4:5)
table(confounders_imp$maternal_edu)

saveRDS(confounders_imp, "confounders_imp_cbclsyndrome_meduRelevel_aromaonly_notwinsib.rds")

c <- readRDS("confounders_imp_cbclsyndrome_meduRelevel_aromaonly_notwinsib.rds")

names(c)
