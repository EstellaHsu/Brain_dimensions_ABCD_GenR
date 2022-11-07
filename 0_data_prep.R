############## ABCD quality control % selection of participants

setwd("/Users/estella/Desktop/ABCD_download")

id <- readRDS("id_has_conMat.rds")
length(id)
#########################################
#demographic info
#########################################

participants <- read.table(file = 'participants.tsv', sep = '\t', header = TRUE)

demo_info <- participants[, c("participant_id","session_id","scanner_manufacturer",
                             "sex","race_ethnicity","age","parental_education")]
demo_info <- demo_info %>% filter(as.character(session_id) == "ses-baselineYear1Arm1")
dim(demo_info)
# there are duplicated ids, but have exact the same information
demo_info <- demo_info[!duplicated(demo_info$participant_id),]
#########################################
#Site info
#########################################

site <- read.delim(file = 'abcd_lt01.txt')
site$participant_id <- sapply(as.character(site$src_subject_id), 
                   function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
site1 <- site %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
site_info <- site1[, c("participant_id","site_id_l")]
names(site_info)[2] <- "site"


#########################################
#genetic: twins & siblings
#########################################

twin <- read.delim("acspsw03.txt")
names(twin)

twin$participant_id <- sapply(as.character(twin$src_subject_id), function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
twin1 <- twin %>% filter(eventname == "baseline_year_1_arm_1")
twin_info <- twin1[,c("participant_id","rel_family_id")]

#########################################
#brain data QC
#########################################

setwd("/Users/estella/Desktop/ABCD_download/manual_qc")
qc_rs <- read.delim("abcd_imgincl01.txt")
names(qc_rs)

qc_rs$participant_id <- sapply(as.character(qc_rs$src_subject_id), 
                        function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
qc1 <- qc_rs %>% filter(eventname == "baseline_year_1_arm_1")

rsfmri_qc <- qc1[,c("participant_id","imgincl_rsfmri_include")]

#########################################
#brain data incidental findings
#########################################

setwd("/Users/estella/Desktop/ABCD_download")
qc_findings <- read.delim("abcd_mrfindings02.txt")
names(qc_findings)

qc_findings$participant_id <- sapply(as.character(qc_findings$src_subject_id), 
                          function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
qc2 <- qc_findings %>% filter(eventname == "baseline_year_1_arm_1")

rsfmri_qc2 <- qc2[,c("participant_id","mrif_score")]


#########################################
#merge all the data
#########################################

#put all data frames into list
df_list <- list(rsfmri_qc, rsfmri_qc2, twin_info, site_info, demo_info)

#merge all data frames in list
all_info <- df_list %>% reduce(full_join, by='participant_id')
dim(all_info)

# 1. id has connectivity matrices
all_info <- all_info %>% filter(participant_id %in% id)

# 2. inclusion criteria
all_info1 <- all_info %>% filter(imgincl_rsfmri_include == 1 & mrif_score < 3)
dim(all_info1)

# 2. remove twin/siblings
all_info2 <- all_info1[!duplicated(all_info1$rel_family_id), ]
dim(all_info2)


# Note: different inclusion criteria has an influence: should start from id has conMat, 
# Otherwise will have less participants (more twin/siblings were excluded)
# the correct sample size should be 6601.

########################################
#merge CBCL
#########################################

setwd("/Users/estella/Desktop/ABCD_download/ABCD_CBCL")
cbcl <- read.delim("abcd_cbcls01.txt")
cbcl$participant_id <- sapply(as.character(cbcl$src_subject_id),function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

cbcl1 <- cbcl %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
names(cbcl1)
######## extract syndrome scale 
cbcl_syn <- cbcl1[, c("participant_id","cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                      "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                      "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                      "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r")]
dim(cbcl_syn)


########################################
#merge all the data
#########################################
all_cbcl <- merge(all_info2, cbcl_syn,by="participant_id")
#str(all_cbcl)
dim(all_cbcl)

all_cbcl1 <- all_cbcl %>% mutate_at(c("age","cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                                      "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                                      "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                                      "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r"),as.numeric) 


setwd("/Users/estella/Desktop/ABCD_download")
id <- readRDS("id_has_conMat.rds")
length(id)
all_cbcl2 <- all_cbcl1[all_cbcl1$participant_id %in% id, ]
dim(all_cbcl2)


colSums(is.na(all_cbcl2))
all_cbcl2[all_cbcl2 == 777] <- NA
all_cbcl2[all_cbcl2 == 888] <- NA
all_cbcl2[all_cbcl2 == 999] <- NA
all_cbcl2[all_cbcl2 == ""] <- NA

colSums(is.na(all_cbcl2))

# remove data with NAS in relevant variables
idx_noNA <- which(rowSums(is.na(all_cbcl2)) == 0)

all_final_noNA <- all_cbcl2[idx_noNA, ]

# change the data type
str(all_final_noNA)

all_final_noNA <- all_final_noNA %>% mutate_at(c("sex","site","race_ethnicity",
                                                 "parental_education","scanner_manufacturer"), as.factor)
str(all_final_noNA)
colSums(is.na(all_final_noNA))

# adjust the level of parental education
levels(all_final_noNA$parental_education) <- list("low" = 1:12,"medium" = 13:17,"high" = 18:21)
str(all_final_noNA)
dim(all_final_noNA)
setwd("/Users/estella/Desktop/ABCD_download/data/newdata_ReleaseQC")
saveRDS(all_final_noNA,"all_final_noNA_incidental_6529.rds")


table(all_final_noNA$site)

dim(all_final_noNA)
########################################
#train test split: site
#########################################
# first split the site 10 times
sites <- unique(all_final_noNA$site)
# there are 22 sites, so 22*0.8=17.6 for the training sites
train_sites <- lapply(1:10, function(i) sample(sites, size = 18, replace = FALSE))

# then split of the ids
all_final_train <- lapply(1:10, function(i) all_final_noNA %>% filter(site %in% train_sites[[i]]))
all_final_test <- lapply(1:10, function(i) all_final_noNA %>% filter(!site %in% train_sites[[i]]))

lapply(all_final_train, nrow)
class(all_final_train[[1]])
saveRDS(all_final_train,"all_final_train.rds")
saveRDS(all_final_test,"all_final_test.rds")

########################################
#ses info: one split
#########################################
setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/retraintest")
train_ses <- readRDS("all_final_train.rds")
test_ses <- readRDS("all_final_test.rds")

summary(train_ses[[1]])
table(train_ses[[1]]$race_ethnicity)/nrow(train_ses[[1]])
table(train_ses[[1]]$parental_education)/nrow(train_ses[[1]])
119.1/12
summary(test_ses[[1]])
table(test_ses[[1]]$race_ethnicity)/nrow(test_ses[[1]])
table(test_ses[[1]]$parental_education)/nrow(test_ses[[1]])

sd(test_ses[[1]]$cbcl_scr_syn_totprob_r)
sd(test_ses[[1]]$age)
########################################
#train test split: site * siemens only
#########################################
all_final_siemens <- all_final_noNA %>% filter(scanner_manufacturer=="SIEMENS")
table(all_final_siemens$site)
# not all sites have siemens
sites_siemens <- sites[table(all_final_siemens$site) > 10]
length(sites_siemens)
# there are 13 sites, so 10 for the training sites
train_sites_siemens <- lapply(1:10, function(i) sample(sites_siemens, size = 10, replace = FALSE))

# then split of the ids
all_final_train_siemens <- lapply(1:10, function(i) all_final_siemens %>% filter(site %in% train_sites_siemens[[i]]))
all_final_test_siemens <- lapply(1:10, function(i) all_final_siemens %>% filter(!site %in% train_sites_siemens[[i]]))

lapply(all_final_train_siemens, nrow)

saveRDS(all_final_train_siemens,"all_final_train_siemens.rds")
saveRDS(all_final_test_siemens,"all_final_test_siemens.rds")















# read the id has the conMat
id <- readRDS("id_has_conMat.rds")
length(id)
all_info0 <- all_info[all_info$participant_id %in% id, ]
#qc1 <- qc[qc$participant_id %in% id, ]
#qc2 <- qc1 %>% filter(as.character(fmri_postqc_qc) == 0)
#dim(qc2)  

# only 1076 has the qc values
#qc2[order(qc2$participant_id), ]
#qc2$participant_id[duplicated(qc2$participant_id)]
#names(qc2)[2] <- "age"
#qc2[qc2$participant_id == "sub-NDARINVP0T27M9X", ]

#id_bad <- qc1$participant_id[qc1$fmri_postqc_qc == 0]
#length(unique(id_bad))
#all_info1 <- all_info0[!all_info0$participant_id %in% id_bad, ]
#dim(all_info1)
#all_info2 <- all_info1[!duplicated(all_info1$participant_id),]
#dim(all_info2)

#all_demo_qc <- all_info2 
#dim(all_demo_qc)

# merge the info from genes
names(gene1)
all_demo_qc_gene <- merge(all_demo_qc, gene1[,c(34,11,14,22)], by="participant_id")
dim(all_demo_qc_gene)

names(all_demo_qc_gen)
#all_demo_qc_gene[order(all_demo_qc_gene$rel_family_id), c(1,8,10, 11,12 )]
all_demo_qc_gen <- all_demo_qc_gene[!duplicated(all_demo_qc_gene$rel_family_id), ]
dim(all_demo_qc_gen)

######### release 4.0: QC criteria

setwd("/Users/estella/Desktop/ABCD_download/manual_qc")
qc_rs <- read.delim("abcd_imgincl01.txt")
head(qc_rs)
qc_rs <- qc_rs[-1, ]

qc_rs$participant_id <- sapply(as.character(qc_rs$src_subject_id), 
                            function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
qc2 <- qc_rs[qc_rs$participant_id %in% id, ]
dim(qc2)

qc3 <- qc2 %>% filter(eventname == "baseline_year_1_arm_1")
dim(qc3)
names(qc3)
id_bad1 <- qc3$participant_id[qc3$imgincl_rsfmri_include == 0]
length(unique(id_bad1))

all_info0 <- all_info[all_info$participant_id %in% id, ]
all_info1 <- all_info0[!all_info0$participant_id %in% id_bad1, ]
dim(all_info1)
all_info1[all_info1$participant_id=="sub-NDARINV007W6H7B",]
all_info2 <- all_info1[!duplicated(all_info1$participant_id),]
dim(all_info2)
all_demo_qc <- all_info2 
dim(all_demo_qc)

# merge the info from genes
names(gene1)
all_demo_qc_gene <- merge(all_demo_qc, gene1[,c(34,11,14,22)], by="participant_id")
dim(all_demo_qc_gene)

names(all_demo_qc_gen)
#all_demo_qc_gene[order(all_demo_qc_gene$rel_family_id), c(1,8,10, 11,12 )]
all_demo_qc_gen <- all_demo_qc_gene[!duplicated(all_demo_qc_gene$rel_family_id), ]
dim(all_demo_qc_gen)



#apply(all_demo_qc_gen[, -1], 2, table)
#a <- all_demo_qc_gen[all_demo_qc_gen$rel_relationship == 1,10:12]
#a[order(a$rel_family_id),]

#saveRDS(all_demo_qc_gen, "all_demo_qc_removeTWINSIB.rds")

# add the gene info
#gene <- twins_sib_abcd %>% filter(participant_id %in% all_demo_qc$participant_id)
#dim(gene)
#gene[order(gene$participant_id), ]
#gene1 <- gene %>% filter(!rel_relationship == "")
#gene1[order(gene1$rel_family_id), ]

#length(id)
#all_demo_qc[order(all_demo_qc$participant_id), ]

#idd <- all_info1$participant_id
#idd[duplicated(idd)]
#all_info1[all_info1$participant_id == "sub-NDARINV1R5PJCK4", ]

######## remove missingness of covariates
colSums(is.na(all_demo_qc_gen)) # nomissingness

######## remove missingness of CBCL: item level
setwd("/Users/estella/Desktop/ABCD_download/ABCD_CBCL")
cbcl <- read.delim("abcd_cbcl01.txt")
cbcl <- cbcl[-1, ]
names(cbcl)

cbcl$participant_id <- sapply(as.character(cbcl$src_subject_id), 
                              function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

cbcl1 <- cbcl %>% filter(participant_id %in% all_demo_qc_gen$participant_id & as.character(eventname) == "baseline_year_1_arm_1")
dim(cbcl1)
names(cbcl1)
dim(all_demo_qc_gen)
cbcl2 <- cbcl1[, c(132, 10:128)]
# merge all the info
all_demo_qc_gen_cbcl <- merge(all_demo_qc_gen,cbcl2,by="participant_id")
dim(all_demo_qc_gen_cbcl)

#data.frame(A=id[id%in% all_demo_qc_gen_cbcl$participant_id], B=all_demo_qc_gen$participant_id)
#ids <- id[id%in% all_demo_qc_gen_cbcl$participant_id]
identical(as.vector(ids), as.vector(all_demo_qc_gen$participant_id))

all_demo_qc_gen_cbcl[all_demo_qc_gen_cbcl == 777] <- NA
all_demo_qc_gen_cbcl[all_demo_qc_gen_cbcl == 888] <- NA
all_demo_qc_gen_cbcl[all_demo_qc_gen_cbcl == 999] <- NA
names(all_demo_qc_gen_cbcl)


colSums(is.na(all_demo_qc_gen_cbcl))
dim(all_demo_qc_gen_cbcl)
length(which(rowSums(is.na(all_demo_qc_gen_cbcl[,c("sex","race_ethnicity","parental_education","rel_family_id")]))!=0))
# remove the kids with covariate missing
all_demo_qc_gen_cbcl_item <- all_demo_qc_gen_cbcl[rowSums(is.na(all_demo_qc_gen_cbcl[,c("sex","race_ethnicity","parental_education","rel_family_id")])) == 0, ]
dim(all_demo_qc_gen_cbcl_item)
colSums(is.na(all_demo_qc_gen_cbcl_item))



all_demo_qc_gen_cbcl_item$sex <- as.factor(all_demo_qc_gen_cbcl_item$sex)
all_demo_qc_gen_cbcl_item$race_ethnicity <- as.factor(all_demo_qc_gen_cbcl_item$race_ethnicity)
all_demo_qc_gen_cbcl_item$parental_education <- as.factor(all_demo_qc_gen_cbcl_item$parental_education)
levels(all_demo_qc_gen_cbcl_item$parental_education) <- list("low" = 1:12,
                                                               "medium" = 13:17,
                                                               "high" = 18:21)
str(all_demo_qc_gen_cbcl_item)
sum(table(all_demo_qc_gen_cbcl_item$site))

names(all_demo_qc_gen_cbcl_item)
all_demo_qc_gen_cbcl_item[, 13:131] <- lapply(all_demo_qc_gen_cbcl_item[, 13:131], function(x) as.numeric(as.character(x)))
str(all_demo_qc_gen_cbcl_item)

setwd("/Users/estella/Desktop/ABCD_download/data/newdata_withManualQC/item_level")
saveRDS(all_demo_qc_gen_cbcl_item, "all_demo_qc_gene_cbclitem.rds")


################### select the cbcl items that I used
names(all_demo_qc_gen_cbcl_new)
cbcl2 <- all_demo_qc_gen_cbcl_new[, 13:131]
dim(cbcl2)
cbcl_abcd <- cbcl2[, names(cbcl2) %in% c("cbcl_q105_p", "cbcl_q101_p", "cbcl_q106_p", "cbcl_q18_p", "cbcl_q99_p",
                                         "cbcl_q72_p","cbcl_q73_p", "cbcl_q97_p", "cbcl_q02_p", "cbcl_q110_p")]


v <- apply(cbcl2, 2, function(x){var(x, na.rm = TRUE)})
v[order(v)]
dim(cbcl_abcd)
colSums(is.na(cbcl_abcd))
sum(rowSums(cbcl_abcd == "") != 0)

apply(all_demo_qc_gen_cbcl_new[,-c(1,2,10)],2,table)
hist(all_demo_qc_gen_cbcl_new$age)
sum(all_demo_qc_gen_cbcl_new$rel_family_id == "")

cbcl_abcd$participant_id <- all_demo_qc_gen_cbcl_new$participant_id

saveRDS(cbcl_abcd, "cbcl_abcd_0102.rds")

str(cbcl_abcd)
cbcl_abcd[1921,]
dim(all_demo_qc_gen_cbcl_new)


dim(all_demo_qc_gen)


#################################
########### CBCL sum scores #####
#################################
setwd("/Users/estella/Desktop/ABCD_download/ABCD_CBCL")
cbcl <- read.delim("abcd_cbcls01.txt")
cbcl <- cbcl[-1, ]
names(cbcl)
cbcl$participant_id <- sapply(as.character(cbcl$src_subject_id), 
                              function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

cbcl1 <- cbcl %>% filter(participant_id %in% all_demo_qc_gen$participant_id & as.character(eventname) == "baseline_year_1_arm_1")
# extract syndrome scale 
cbcl_syn <- cbcl1[, c(92, 10:53)]
dim(cbcl_syn)

names(cbcl_syn)
cbcl_syn_r <- cbcl_syn[, c(seq(2, 42, 4))]
dim(cbcl_syn_r)
cbcl_syn_r[cbcl_syn_r == 999] <- NA
cbcl_syn_r[cbcl_syn_r == 888] <- NA
cbcl_syn_r[cbcl_syn_r == 777] <- NA
cbcl_syn_r[cbcl_syn_r == ""] <- NA
# the original data is character, now need to be changed to numbers
cbcl_syn_r[] <- lapply(cbcl_syn_r, function(x) {as.numeric(as.character(x))})
which(colSums(is.na(cbcl_syn_r)) != 0)
cbcl_syn_rawscore <- cbcl_syn_r
cbcl_syn_rawscore$participant_id <- cbcl_syn$participant_id
names(cbcl_syn_rawscore)
dim(cbcl_syn_rawscore)
setwd("/Users/estella/Desktop/ABCD_download/data/syndrome_scale/")
saveRDS(cbcl_syn_rawscore, "cbcl_syn_rawscore_abcd.rds")

all_demo_qc_gen_cbclRawscore <- merge(all_demo_qc_gen,cbcl_syn_rawscore,by="participant_id")
dim(all_demo_qc_gen_cbclRawscore)
all_demo_qc_gen_cbclRawscore[all_demo_qc_gen_cbclRawscore == 999] <- NA
all_demo_qc_gen_cbclRawscore[all_demo_qc_gen_cbclRawscore == 888] <- NA
all_demo_qc_gen_cbclRawscore[all_demo_qc_gen_cbclRawscore == 777] <- NA
colSums(is.na(all_demo_qc_gen_cbclRawscore))
names(all_demo_qc_gen_cbclRawscore)
summary(all_demo_qc_gen_cbclRawscore)
# remove NAs
all_demo_qc_gen_cbclRawscore_final <- all_demo_qc_gen_cbclRawscore[rowSums(is.na(all_demo_qc_gen_cbclRawscore)) == 0, ]
dim(all_demo_qc_gen_cbclRawscore_final)
colSums(is.na(all_demo_qc_gen_cbclRawscore_final))

# adjust the structure of the data
all_demo_qc_gen_cbclRawscore_final$sex <- as.factor(all_demo_qc_gen_cbclRawscore_final$sex)
all_demo_qc_gen_cbclRawscore_final$race_ethnicity <- as.factor(all_demo_qc_gen_cbclRawscore_final$race_ethnicity)
all_demo_qc_gen_cbclRawscore_final$parental_education <- as.factor(all_demo_qc_gen_cbclRawscore_final$parental_education)
levels(all_demo_qc_gen_cbclRawscore_final$parental_education) <- list("low" = 1:12,
                                                             "medium" = 13:17,
                                                             "high" = 18:21)
str(all_demo_qc_gen_cbclRawscore_final)
sum(table(all_demo_qc_gen_cbclRawscore_final$sex))

setwd("/Users/estella/Desktop/ABCD_download/data/newdata_withManualQC/syndrome_scale")
saveRDS(all_demo_qc_gen_cbclRawscore_final, "all_demo_qc_gene_cbclsyndromeRawscore_newQC.rds")


cbcl_item_genr <- cbcl_ite_genr[, !names(cbcl_ite_genr) %in% c("cbcl105_9m", "cbcl101_9m", "cbcl106_9m", "cbcl18_9m",
                                                               "cbcl72_9m","cbcl73_9m", "cbcl97_9m",  "cbcl110_9m")]

cbcl_abcd <- cbcl2[, names(cbcl2) %in% c("cbcl_q105_p", "cbcl_q101_p", "cbcl_q106_p", "cbcl_q18_p", "cbcl_q99_p",
                                         "cbcl_q72_p","cbcl_q73_p", "cbcl_q97_p", "cbcl_q02_p", "cbcl_q110_p")]










