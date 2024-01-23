
############################################################################
############## ABCD quality control % selection of participants ############
############################################################################

# the ids that have the connecitvity matrices
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

qc_rs <- read.delim("abcd_imgincl01.txt")

qc_rs$participant_id <- sapply(as.character(qc_rs$src_subject_id), 
                        function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
qc1 <- qc_rs %>% filter(eventname == "baseline_year_1_arm_1")

rsfmri_qc <- qc1[,c("participant_id","imgincl_rsfmri_include")]

#########################################
#brain data QC 2: framewise displacement
#########################################
                               
fd_abcd <- read.delim("abcd_betnet02.txt")
fd_abcd$participant_id <- sapply(as.character(fd_abcd$src_subject_id), 
                                 function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
fd1 <- fd_abcd %>% filter(as.character(eventname) == "baseline_year_1_arm_1")

fd_info <- fd1[, c("participant_id","rsfmri_c_ngd_meanmotion")]
fd_info$rsfmri_c_ngd_meanmotion <- as.numeric(fd_info$rsfmri_c_ngd_meanmotion)

                             
#########################################
#brain data incidental findings
#########################################

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
df_list <- list(rsfmri_qc, rsfmri_qc2, twin_info, site_info, demo_info, fd_info)

#merge all data frames in list
all_info <- df_list %>% reduce(full_join, by='participant_id')
dim(all_info)

# 1. id has connectivity matrices
all_info <- all_info %>% filter(participant_id %in% id)
dim(all_info)

# 2. inclusion criteria
all_info1 <- all_info %>% filter(imgincl_rsfmri_include == 1 & mrif_score < 3 & rsfmri_c_ngd_meanmotion < 0.25)
dim(all_info1)

# 3. remove twin/siblings
all_info2 <- all_info1[!duplicated(all_info1$rel_family_id), ]
dim(all_info2)

                                                                    
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

all_cbcl1[all_cbcl1 == 777] <- NA
all_cbcl1[all_cbcl1 == 888] <- NA
all_cbcl1[all_cbcl1 == 999] <- NA
all_cbcl1[all_cbcl1 == ""] <- NA

colSums(is.na(all_cbcl1))

# remove data with NAS in relevant variables
idx_noNA <- which(rowSums(is.na(all_cbcl2)) == 0)

all_final_noNA <- all_cbcl2[idx_noNA, ]

# change the data type
str(all_final_noNA)

all_final_noNA <- all_final_noNA %>% mutate_at(c("sex","site","race_ethnicity",
                                                 "parental_education","scanner_manufacturer"), as.factor)
str(all_final_noNA)


# adjust the level of parental education
levels(all_final_noNA$parental_education) <- list("low" = 1:12,"medium" = 13:17,"high" = 18:21)
str(all_final_noNA)
dim(all_final_noNA)

saveRDS(all_final_noNA,"all_final_noNA_incidental_6529.rds")


########################################
#train test split: site
#########################################
# first split the site 30 times
sites <- unique(all_final_noNA$site)
# there are 18 sites for the training sites
train_sites <- lapply(1:30, function(i) sample(sites, size = 18, replace = FALSE))

# then split of the ids
all_final_train <- lapply(1:30, function(i) all_final_noNA %>% filter(site %in% train_sites[[i]]))
all_final_test <- lapply(1:30, function(i) all_final_noNA %>% filter(!site %in% train_sites[[i]]))
                     
saveRDS(all_final_train,"all_final_train.rds")
saveRDS(all_final_test,"all_final_test.rds")


