#################################################################
###################### ABCD: code for figures ###################
#################################################################

# These are the code for the figures used in the manuscript

########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'fmsb', 'NbClust','stats','dendextend','cluster',
              'fpc', 'e1071', 'plotly', 'MASS', 'ggradar2','tidyr','introdatavizyes','Hmisc',
              'pheatmap','RColorBrewer','circlize')


lapply(packages, require, character.only = TRUE)
###############################################################
############## Visualization of correlations across 10 splits 
###############################################################
####### reorder the canonical variates
rs_train_test_abcd <- readRDS("rs_train_test_abcd.rds")
# the reference split: choose one split as the reference split
load_std <- res.abcd.ref$v # the loadings for CBCL

new_cv_reorder  <- lapply(seq_along(rs_train_test_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:6) {
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:6, function(z) {cor(load_std[, i], rs_train_test_abcd[[j]]$abcd.train$v[, z])})
    if(max(abs(cor_cv)) < 0.5){ # sometimes the correlation is too low
      idx  <- i
    } else {
      idx <- which(abs(cor_cv) == max(abs(cor_cv)))
    }
    
    if(!idx %in% idx_v) {
      idx_v <- c(idx_v, idx) # sometimes it will overlap
    } else {
      idx_v <- c(idx_v, i) 
    }
  }
  idx_v[duplicated(idx_v)] <- c(1:6)[!c(1:6) %in% idx_v]
  
  loadMat <- rs_train_test_abcd[[j]]$abcd.train$v[, idx_v]
  return(list(loadMat=loadMat, idx_reorder=idx_v))
  
})


###### reorder the correlations
traincor <- lapply(rs_train_test_abcd, function(x) {x$abcd.train$cors})
traincor <- do.call(rbind,traincor)
traincor_reorder <- lapply(1:10, function(i) {traincor[i, new_cv_reorder[[i]]$idx_reorder]})
traincor_reorder <- do.call(rbind,traincor_reorder)
rownames(traincor_reorder) <- paste0("split",1:10)
colnames(traincor_reorder) <- paste0("CV",1:6)
meantrain <- colMeans(traincor_reorder)
sdtrain <- apply(traincor_reorder, 2, sd)


testcor <- lapply(rs_train_test_abcd, function(x) {x$abcd.test})
testcor <- do.call(rbind,testcor)
testcor_reorder <- lapply(1:10, function(i) {testcor[i, new_cv_reorder[[i]]$idx_reorder]})
testcor_reorder <- do.call(rbind,testcor_reorder)
rownames(testcor_reorder) <- paste0("split",1:10)
colnames(testcor_reorder) <- paste0("CV",1:6)
meantest <- colMeans(testcor_reorder)
sdtest <- apply(testcor_reorder, 2, sd)

fig_10splits_syn <- data.frame(Train_Test = rep(c("Training set","Test set"),each=3), meancor = c(meantrain[1:3], meantest[1:3]),
                               sdcor = c(sdtrain[1:3], sdtest[1:3]), CV=rep(paste0("CV",1:3),2))


p_cor <- ggplot(fig_10splits_syn, aes(x=CV, y=meancor, color=Train_Test, group=Train_Test)) + 
  scale_color_manual(name="Training Test Set", values = c("#fbb172", "#E64659")) + 
  geom_errorbar(aes(ymin=meancor-sdcor, ymax=meancor+sdcor), width=.2, size=2) +
  geom_point(size=5) + ylim(c(0, 0.25)) + 
  theme_bw() + labs(y = "canonical correlations", x = "canonical variates", size=5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text = element_text(size=13))
p_cor

ggsave("train_test_correlations.pdf", width = 8.5, height = 6)



###############################################################
############## Radar plot of CBCL syndrome loadings  ##########
###############################################################

loadings_cv <- lapply(new_cv_reorder, function(x) {x$loadMat})

# calculate the mean of across the 10 splits
temp <- sapply(1:6, function(i) {
  rowMeans(sapply(seq_along(loadings_cv), function(j) {abs(loadings_cv[[j]][, i])}))
})

df_temp <- as.data.frame(temp)

##################################
######### circular radar plot
##################################
df_temp1 <- t(df_temp)
df_temp2 <- rbind(rep(1,8), df_temp1)
colnames(df_temp2) <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule breaking","Aggression")

df_temp3 <- as.data.frame(cbind(group=c("1","CV1", "CV2", "CV3", "CV4", "CV5", "CV6"),df_temp2))
df_temp3[,2:9] <- lapply(df_temp3[,2:9], function(x){as.numeric(as.character(x))})

colors_border_6 <- c("#fcfcfb","#A0D568","#F7EA48","#4FC1E8","#F08080","#CC99C9","#FFCE54")
colors_in_6 <- c("#fcfcfb","#A0D568","#F7EA48","#4FC1E8","#F08080","#CC99C9","#FFCE54")

p2 <- ggradar2(df_temp3, base.size=5, webtype = 'lux', grid.min = 0,
         grid.max = 1,label.gridline.mid = TRUE, group.colours = colors_border_6,
         group.fill.colours = colors_in_6,label.centre.y=FALSE,
         gridline.mid.colour="grey", grid.label.size = 0,
         gridline.max.linetype = "solid",polygonfill.transparency=0.5,
         background.circle.transparency=0.1,axis.label.size=5.5,
         group.line.width=1.5,group.point.size=3,plot.legend = TRUE) 
p2

ggsave("radarfigure_train_abcd.pdf", width = 12, height = 6)

dev.off()


########################################################
############## Covariance explained  ###################
########################################################

# This is an example from one of the 10 train-test splits
# example for the training set
# code for the test set is the same, just need to change the "brain_train" to "brain_test", "cbcl_train" to "cbcl_test"
vardf <- VarianceExplain(brain_train, cbcl_train, res.abcd.example, 8) 
colnames(vardf) <- c("Canonical_Variates", "Covariance_Explained")
vardf$`Canonical_Variates` <- paste0("CV",1:8)
vardf$Canonical_Variates <- as.factor(vardf$Canonical_Variates)

ggplot(vardf, aes(x=Canonical_Variates, y=Covariance_Explained, group=Canonical_Variates))+
  geom_point(aes(shape=Canonical_Variates,color=Canonical_Variates), size = 7, alpha = 0.9) + 
  theme_bw() + theme(legend.position="none",
                     axis.text.x = element_text(size = 12)) + 
  scale_shape_manual(values=c(17,16,15,18,20,20,8,8)) +
  scale_color_manual(values=c("#A0D568","#F7EA48","#4FC1E8","#f08080","#cc99c9","#ffce54","grey","grey")) 

ggsave("covariance_explained_trainset.pdf",width=8,height=4)



########################################################
############## permutation test visualization  #########
########################################################

# This is an example from one of the 10 train-test splits
# example for the test set

cor.abcd.test.example <- test_project_weights(brain_test, cbcl_test, res.abcd.example, 8)
perm_abcdtest.example <- permutation_test_testset(cbcl_test, brain_test,nperm=1999, 
                                           res.abcd.example,cor.abcd.test.example)


# fdr correction for the p values
p.fdr <- p.adjust(perm_abcdtest.example$pval.perm, method="fdr")

# example for the first canonical correlation
# this is the null distribution of the permutation test
df.cor.1 <- data.frame(cor = perm_abcdtest.example$cor.perm[, 1]) 

pdf("permutation_test_cv1.pdf", width=5, height=3)
ggplot(df.cor.1, aes(x=cor)) + 
  geom_histogram(binwidth=0.01, fill = "#4FC1E8", color = "black", alpha = 0.8)+
  geom_vline(xintercept=p.fdr[1], linetype="dashed", color = "red", size = 1) +
  ggtitle("Permutation correlations of the 1st canonical variate") + 
  labs(x = "canonical correlations") + 
  xlim(c(-0.2, 0.2)) + annotate(geom="text", x=0.15, y=180, label="P(FDR) < 0.01", 
                                fontface = 'italic', size = 5, color="black") + 
  theme_bw()
dev.off()



#######################################################################
############################ Brain figures  ###########################
#######################################################################


########################################
######## Read the data 
########################################
# This is calculated using weighted PCA in the whole ABCD sample

pca_whole <- readRDS("pca.weighted.whole.rds")
str(pca_whole)
brain_whole <- pca_whole$brain_train_reduced[, 1:100] # the first 100 brain PCs

setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC")
all_whole <- readRDS("all_final_noNA_incidental_6529.rds")
cbcl_whole <- all_whole[, c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                        "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                        "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")]


########################################
######## Calculate brain scores
########################################
# brain mean is the average brain loadings across 10 splits
brain_cvscores <- scale(brain_whole) %*% brain_mean

# Then calculate the correlations between the original brain data and the brain canonical scores
# This is used as the indicator of the importance of each brain feature
# "feature_abcd" is the original vectorized brain features
feature_abcd <- readRDS("brain_whole_residual.rds")
feature_abcd <- as.data.frame(feature_abcd)
cor_brain_cv <- lapply(1:ncol(feature_abcd), function(i) {cor(feature_abcd[, i], brain_cvscores)})
cor_brain_cv <- do.call(rbind,cor_brain_cv)

corBrCv1 <- cor_brain_cv[, 1]
corBrCv2 <- cor_brain_cv[, 2]
corBrCv3 <- cor_brain_cv[, 3]

########################################
######## read the brain labels
########################################
gordon <- read.table("Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.CLUT.txt", sep = "\t")
labels <- gordon[-seq(2, 704, 2), ]
length(labels)
labels <- gsub('^[LR]_', "", labels)
labels_group <- gsub("\\_.*", "",labels)

# labels for the 352 * 352 matrices
labels_group_sub <- c(labels_group[1:333],rep("Subcortical", 19))
# labels for the parcel-based matrices
parcels_sub <- unique(labels_group_sub)

Nparcels <- sapply(unique(labels_group_sub), function(x) {sum(labels_group_sub == x)})
parcels <- unique(labels_group)


######################################
####### create the heatmap
######################################
This is an example 

# build an empty brain matrix
brainMat <- matrix(NA,352,352)
# fill in the brain loadings for the first canonical variate (others are the same)
brainMat[upper.tri(brainMat)] <- corBrCv1
brainMat[lower.tri(brainMat, diag = F)] <- t(brainMat)[lower.tri(brainMat)]
colnames(brainMat) <- labels_group_sub
rownames(brainMat) <- labels_group_sub

# build an empty matrix for the parcels-based matrix
group_connect <- matrix(0,length(parcels_sub),length(parcels_sub))
colnames(group_connect) <- parcels_sub
rownames(group_connect) <- parcels_sub

m_connect <- brainMat
# calculate within and between modules connectivity: equations based on Xia et al. (2018)

fig_df <- lapply(seq_along(parcels_sub), function(i){
  
  idxr <- which(rownames(m_connect) == parcels_sub[i])
  M <- sum(rownames(m_connect) == parcels_sub[i])
  within_con <- m_connect[idxr, idxr]
  within_load <- 2*sum(within_con, na.rm = TRUE)/(M * (M-1))
  
  if(i < 14){ 
    between_con <- lapply((i+1):length(parcels_sub), function(j) {
      idxc <-  which(colnames(m_connect) == parcels_sub[j])
      between <- m_connect[idxr, idxc]
      M <- sum(rownames(m_connect) == parcels_sub[i])
      N <- sum(colnames(m_connect) == parcels_sub[j])
      list(M = M, N = N, between_con = between, from = parcels_sub[i], to = parcels_sub[j])})
  } else {between_con <- NULL}
  
  between_load <- data.frame()
  for (j in seq_along(between_con)){
    M <- unlist(between_con[[j]][1])
    N <- unlist(between_con[[j]][2])
    conM <- between_con[[j]][3]
    bet_load <- sum(unlist(conM))/(M * N)
    from <- unlist(between_con[[j]][4])
    to <- unlist(between_con[[j]][5])
    temp <- data.frame(loadings = bet_load, from = from, to = to)
    between_load <- rbind(between_load, temp)
  }
  
  df1 <- data.frame(loadings = within_load, from = parcels_sub[i], to = parcels_sub[i])
  df2 <- rbind(df1, between_load)
  list(within_load = within_load, between_load = between_load, fig = df2)
  
})


# extract between and within loadings
within <- do.call(rbind, lapply(fig_df, function(x) {x$within_load}))
between <- do.call(rbind, lapply(fig_df, function(x) {x$between_load}))

# calculate the z scores
group_connect[lower.tri(group_connect)] <- scale(between$loadings)
group_fig <- t(group_connect)
group_fig[lower.tri(group_fig)] <- scale(between$loadings)
diag(group_fig) <- scale(within)

breaksList <- seq(-2, 2, by = 0.1)
length(breaksList)
cols <- c(colorRampPalette(c("#4575b4","#abd9e9"))(12), rep("white", 18), colorRampPalette(c("#eff121","#e93131"))(11))


pdf("brain_heatmap_cv1_average10.pdf", width=5.5, height=5)

pheatmap(group_fig, treeheight_row = 0, treeheight_col = 0, 
         cluster_rows = F, cluster_cols = F,
         color = cols,
         breaks = breaksList) # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList))

dev.off()


#########################################################
################### connectogram 
#########################################################

fig_temp <- do.call(rbind, lapply(fig_df, function(x) {x$fig}))
fig_connecto <- as.data.frame(fig_temp)

fig_con <- data.frame(from=fig_connecto$from, to=fig_connecto$to,
                      value=scale(fig_connecto$loadings))
higher <- quantile(fig_con$value, probs = .975)
lower <- quantile(fig_con$value, probs = .025)

# just present the top 5% important brain networks
# < 97.5% or > 2.5%
fig_con[fig_con < higher & fig_con > lower] <- 0

col <- c("Default"="#d55365", "Salience"="#afe0a9", "VentralAttn"="#fbb172",
         "Auditory"="#93C6F6", "Subcortical"="#83DFEA", "ParietoOccip"= "#77c6a9",
         "FrontoParietal"="#5390bf", "MedialParietal"="#FFE43F","DorsalAttn"="#B097D6",
         "Visual"="#a2335d", "SMhand"="#B7E457", "SMmouth"="#F9F34B","None"="grey")


pdf("brain_connectogram_cv1_average10.pdf", width=8, height=8)
par(cex = 1.0, mar = c(0, 0, 0, 0))
chordDiagram(fig_con, grid.col = col, annotationTrack = c("name","grid"),
             annotationTrackHeight = c(0.01, 0.08))
dev.off()

