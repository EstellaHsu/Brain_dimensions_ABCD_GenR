###############################################################################################
############## ABCD: test the associations between brain CVs and cognition ####################
###############################################################################################

############ read the information of all the included participants
allinfo <- readRDS("all_final_noNA_incidental_6529.rds")
id <- allinfo$participant_id

############# read the NIH tool box data
nih_tool <- read.delim("abcd_tbss01.txt")
names(nih_tool)
dim(nih_tool)
nih_tool$participant_id <- sapply(as.character(wisc$src_subject_id), 
                            function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
# filter the baseline data
nih_tool <- nih_tool %>% filter(eventname == "baseline_year_1_arm_1")

nih_abcd <- nih_tool[nih_tool$participant_id %in% id, c("participant_id", "nihtbx_fluidcomp_agecorrected", 
                                             "nihtbx_cryst_agecorrected",
                                             "nihtbx_totalcomp_agecorrected")]
dim(nih_abcd)

############# read the WISC data
wisc <- read.delim("abcd_ps01.txt")
names(wisc)

wisc$participant_id <- sapply(as.character(wisc$src_subject_id), 
                              function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

wisc <- wisc %>% filter(eventname == "baseline_year_1_arm_1")

wisc_abcd <- wisc[wisc$participant_id %in% id, c("participant_id","pea_wiscv_tss")]
dim(wisc_abcd)

########################################################
############## Calculate canonical scores ##############
########################################################
# calculate the brain CV scores based on the average brain canonical weights across 30 splits
CV_abcd_total <- scale(brain_whole) %*% brain_mean

df_cvs <- data.frame(participant_id=id, cv1=CV_abcd_total[,1], 
                     cv2=CV_abcd_total[,2])


########################################################
############## Merge the data ##########################
########################################################

all0 <- merge(df_cvs, nih_abcd, by="participant_id")
dim(all0)
all1 <- merge(all0, wisc_abcd, by="participant_id")
dim(all1)
all2 <- merge(all1, allinfo[,c("participant_id","sex","site","race_ethnicity","age","parental_education")], by="participant_id")
colSums(is.na(all2))
# no missingness in the data                           
all_noNA <- all2
dim(all_noNA)
all_noNA <- all_noNA %>% mutate(matrix=as.numeric(as.character(all_noNA$pea_wiscv_tss),
                                fluid=as.numeric(as.character(all_noNA$nihtbx_fluidcomp_agecorrected),
                                cryst=as.numeric(as.character(all_noNA$nihtbx_cryst_agecorrected),
                                total=as.numeric(as.character(all_noNA$nihtbx_totalcomp_agecorrected))
str(all_noNA)

# last check to remove participants with missingness
all_noNA <- all_noNA[rowSums(is.na(all_noNA[,c("matrix","fluid","cryst","total")])) == 0, ]
dim(all_noNA)

########################################################
############## linear regressions ######################
########################################################
# small function for regressions
cognition <- function(cv, outcome) {
   model <- formula(paste0("scale(",outcome,") ~ ", "scale(", cv, ")+site+sex+race_ethnicity+age+parental_education"))
   lm_model <- lm(model, all_noNA)
   return(lm_model)
 }
                                                                                                  
matrix_cv1 <- cognition("cv1", "matrix")
matrix_cv2 <- cognition("cv2", "matrix")                                                
# make the html table
tab_model(matrix_cv1,matrix_cv2,matrix_cv3)                                             

fluid_cv1 <- cognition("cv1", "fluid")
fluid_cv2 <- cognition("cv2", "fluid")                                                                                               
tab_model(fluid_cv1,fluid_cv2,fluid_cv3)  
                                                 
cryst_cv1 <- cognition("cv1", "cryst")
cryst_cv2 <- cognition("cv2", "cryst")                                                                                               
tab_model(fluid_cv1,fluid_cv2,fluid_cv3) 
                                                 
total_cv1 <- cognition("cv1", "total")
total_cv2 <- cognition("cv2", "total")                                                                                                
tab_model(fluid_cv1,fluid_cv2,fluid_cv3)                                                   
                                                 
