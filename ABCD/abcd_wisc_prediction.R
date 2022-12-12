########################################################
############## ABCD wisc prediction ####################
########################################################


setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/")

allinfo <- readRDS("all_final_noNA_incidental_6529.rds")
id <- allinfo$participant_id

############# read the wisc data
setwd("/Users/estella/Desktop/ABCD_download/abcd_wisc")
wisc <- read.delim("abcd_tbss01.txt")
names(wisc)
dim(wisc)

wisc$participant_id <- sapply(as.character(wisc$src_subject_id), 
                            function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
wisc <- wisc %>% filter(eventname == "baseline_year_1_arm_1")

wisc_abcd <- wisc[wisc$participant_id %in% id, c("participant_id", "nihtbx_fluidcomp_agecorrected", 
                                             "nihtbx_cryst_agecorrected",
                                             "nihtbx_totalcomp_agecorrected")]
dim(wisc_abcd)
wisc_abcd$participant_id[duplicated(wisc_abcd)]
# there are duplicated scores
wisc_abcd[wisc_abcd$participant_id == "sub-NDARINV9B3TN6RL",]
wisc_abcd1 <- wisc_abcd[!duplicated(wisc_abcd$participant_id), ]

dim(wisc_abcd1)

############# read the pearson data
setwd("/Users/estella/Desktop/ABCD_download/abcd_pearson")
pearson <- read.delim("abcd_ps01.txt")
names(pearson)

pearson$participant_id <- sapply(as.character(pearson$src_subject_id), 
                              function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

pearson <- pearson %>% filter(eventname == "baseline_year_1_arm_1")

pearson_abcd <- pearson[pearson$participant_id %in% id, c("participant_id","pea_wiscv_tss")]
dim(pearson_abcd)

########################################################
############## Calculate canonical scores ##############
########################################################
CV_abcd_total <- scale(brain_whole) %*% brain_mean
cv1_abcd <- CV_abcd_total[,1]
cv2_abcd <- CV_abcd_total[,2]
cv3_abcd <- CV_abcd_total[,3]

df_cvs <- data.frame(participant_id=id, cv1=cv1_abcd, cv2=cv2_abcd, cv3=cv3_abcd)


########################################################
############## Merge the data ##########################
########################################################

all0 <- merge(df_cvs, pearson_abcd, by="participant_id")
dim(all0)
all1 <- merge(all0, wisc_abcd1, by="participant_id")
dim(all1)
all2 <- merge(all1, allinfo[,c("participant_id","sex","site","race_ethnicity","age","parental_education")], by="participant_id")
colSums(is.na(all2))
all_noNA <- all2
dim(all_noNA)
names(all_noNA)[5:8] <- c("matrix","fluid","cryst","total")
str(all_noNA)
all_noNA$matrix <- as.numeric(as.character(all_noNA$matrix))
all_noNA$fluid <- as.numeric(as.character(all_noNA$fluid))
all_noNA$cryst <- as.numeric(as.character(all_noNA$cryst))
all_noNA$total <- as.numeric(as.character(all_noNA$total))

all_noNA <- all_noNA[rowSums(is.na(all_noNA[,5:8])) == 0, ]
dim(all_noNA)
all_noNA <- all_noNA[all_noNA$cryst != max(all_noNA$cryst), ]# there's an outlier in crystal
dim(all_noNA)
########################################################
############## linear regressions ######################
########################################################
###### CV1
names(all_noNA)
matrix <- lm(scale(matrix) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(matrix)
confint(matrix)
fluid <- lm(scale(fluid) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(fluid)
confint(fluid)
cryst <- lm(scale(cryst) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(cryst)
confint(cryst)
total <- lm(scale(total) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(total)
confint(total)


mat <- cor(df, use="complete.obs")
corrplot(mat)

summary(all_noNA)

ggplot(all_noNA,aes(cv1, cryst)) +
  geom_point() +
  geom_smooth(method='lm') 


p.adjust(c(0.10762,0.22264,0.220253,0.0262,
           0.00341,0.012060,0.000997,0.000183,
           0.0148,0.000105,0.000647,0.000001))





