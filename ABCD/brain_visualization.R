#################################################################
############## ABCD whole sample for brain visualization#########
#################################################################


########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'fmsb', 'NbClust','stats','dendextend','cluster',
              'fpc', 'e1071', 'plotly', 'MASS', 'ggradar2','tidyr','introdatavizyes','Hmisc',
              'pheatmap','RColorBrewer','circlize')


lapply(packages, require, character.only = TRUE)
########################################################
############## 1. Read the data ########################
########################################################
setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/retraintest")

pca_whole <- readRDS("pca.weighted.whole.rds")
str(pca_whole)
brain_whole <- pca_whole$brain_train_reduced[, 1:100]

setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC")
all_whole <- readRDS("all_final_noNA_incidental_6529.rds")
cbcl_whole <- all_whole[, c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                        "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                        "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")]

########################################################
############## 2. Run the sCCA ########################
########################################################

cv.resample <- CV_sampling(cbcl_whole, brain_whole, 100, 0.8)
x_pen <- seq(0.1, 1.0, length.out = 10) # brain
y_pen <- seq(0.1, 1.0, length.out = 10)

grid.wholesample <- grid.search.cor.Testset(cv.resample$brain_train, cv.resample$cbcl_train,
                                              cv.resample$brain_test, cv.resample$cbcl_test,
                                              x_pen, y_pen,nsample=100)



res.wholesample <- CCA(x=brain_whole, z=cbcl_whole, penaltyx = 0.5, penaltyz = 0.5, typex="standard", typez="standard",
                         niter = 20, K=8)
res.wholesample
# visualize the loadings:
cbcl_loading <- res.wholesample$v

#cbcl_loading[, 3:5] <- -cbcl_loading[, 3:5] 
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")
corrplot(t(cbcl_loading)[1:6,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 1.6, number.cex=1)


#########################################################
################### Calculate brain scores
#########################################################
brain_cvscores <- scale(brain_whole) %*% brain_mean
saveRDS(brain_cvscores, "brain_cvscores.rds")
cor_brain <- readRDS("cor_brain_cvscores_wholesample_average10.rds")
cor_brain_cv <- do.call(rbind,cor_brain)

range(cor_brain_cv)

corBrCv1 <- cor_brain_cv[, 1]
corBrCv2 <- cor_brain_cv[, 2]
corBrCv3 <- cor_brain_cv[, 3]

#########################################################
################### Heatmap with Parcellation
#########################################################
setwd("~/Desktop/ABCD_download/")
f <- readRDS("feature_name_reference.rds")
gordon <- read.table("Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.CLUT.txt", sep = "\t")
labels <- gordon[-seq(2, 704, 2), ]
length(labels)
labels <- gsub('^[LR]_', "", labels)
labels_group <- gsub("\\_.*", "",labels)

labels_group_sub <- c(labels_group[1:333],rep("Subcortical", 19))
parcels_sub <- unique(labels_group_sub)

Nparcels <- sapply(unique(labels_group_sub), function(x) {sum(labels_group_sub == x)})
parcels <- unique(labels_group)


######################################
###### heatmap
######################################

brainMat <- matrix(NA,352,352)
brainMat[upper.tri(brainMat)] <- corBrCv3
brainMat[lower.tri(brainMat, diag = F)] <- t(brainMat)[lower.tri(brainMat)]
#diag(brainMat) <- max(brainMat)
colnames(brainMat) <- labels_group_sub
rownames(brainMat) <- labels_group_sub

group_connect <- matrix(0,length(parcels_sub),length(parcels_sub))
colnames(group_connect) <- parcels_sub
rownames(group_connect) <- parcels_sub
m_connect <- brainMat

# calculate within and between modules connectivity

fig_df <- lapply(seq_along(parcels_sub), function(i){
  
  idxr <- which(rownames(m_connect) == parcels_sub[i])
  M <- sum(rownames(m_connect) == parcels_sub[i])
  within_con <- m_connect[idxr, idxr]
  within_load <- 2*sum(within_con, na.rm = TRUE)/(M * (M-1))
  
  if(i < 14){ # when i = 14, j will be 15. 
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

fig_temp <- do.call(rbind, lapply(fig_df, function(x) {x$fig}))

fig_connecto <- as.data.frame(fig_temp)


within <- do.call(rbind, lapply(fig_df, function(x) {x$within_load}))
between <- do.call(rbind, lapply(fig_df, function(x) {x$between_load}))

group_connect[lower.tri(group_connect)] <- scale(between$loadings)
group_fig <- t(group_connect)
group_fig[lower.tri(group_fig)] <- scale(between$loadings)
diag(group_fig) <- scale(within)

breaksList <- seq(-2, 2, by = 0.1)
length(breaksList)
cols <- c(colorRampPalette(c("#4575b4","#abd9e9"))(12), rep("white", 18), colorRampPalette(c("#eff121","#e93131"))(11))

setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC/retraintest")

pdf("brain_heatmap_cv1_average10.pdf", width=5.5, height=5)

pheatmap(group_fig, treeheight_row = 0, treeheight_col = 0, 
         cluster_rows = F, cluster_cols = F,
         color = cols,
         breaks = breaksList) # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList))

dev.off()
range(group_fig)


#########################################################
################### connectogram 
#########################################################
library("circlize")
library("tidyverse")

fig_con <- data.frame(from=fig_connecto$from, to=fig_connecto$to,
                      value=scale(fig_connecto$loadings))
quantile(fig_con$value, probs = .95)

fig_con[fig_con < 1.1548 & fig_con > -1.1548] <- 0

col <- c("Default"="#d55365", "Salience"="#afe0a9", "VentralAttn"="#fbb172",
         "Auditory"="#93C6F6", "Subcortical"="#83DFEA", "ParietoOccip"= "#77c6a9",
         "FrontoParietal"="#5390bf", "MedialParietal"="#FFE43F","DorsalAttn"="#B097D6",
         "Visual"="#a2335d", "SMhand"="#B7E457", "SMmouth"="#F9F34B","None"="grey")


pdf("brain_connectogram_cv3_average10.pdf", width=8, height=8)
par(cex = 1.0, mar = c(0, 0, 0, 0))
chordDiagram(fig_con, grid.col = col, annotationTrack = c("name","grid"),
             annotationTrackHeight = c(0.01, 0.08))
dev.off()





















nameMat <- matrix(0, 352,352)
colnames(nameMat) <- labels
rownames(nameMat) <- labels

for (i in 1:351) {
  for (j in (i+1):352) {
    nameMat[i,j] <- paste0(rownames(nameMat)[i], "/", colnames(nameMat)[j])
  }
}

nameMat_group <- matrix(0, 352,352)
colnames(nameMat_group) <- labels_group
rownames(nameMat_group) <- labels_group

for (i in 1:351) {
  for (j in (i+1):352) {
    nameMat_group[i,j] <- paste0(rownames(nameMat_group)[i], "/", colnames(nameMat_group)[j])
  }
}

name_mat <- nameMat[upper.tri(nameMat,diag = FALSE)]
# the vertical version of the name-feature reference for extract the data frame of edges
feature_name_reference <- data.frame(feature=paste0("F_", 1:length(name_mat)), brain_name=name_mat)


name <- lapply(seq(name_mat), function(x) {unlist(strsplit(name_mat[x], "/"))})
from <- unlist(lapply(name, function(x) { return(x[1])}))
to <- unlist(lapply(name, function(x) { return(x[2])}))

# extract the weights
#wcca <- train.pca$u[, 1]
#idx_top10 <- which(order(abs(wcca), decreasing = T) <= 10)
#wpca <- pca.train.weighted$rotation[,1:200]

#new_weightMat <- apply(wpca, 1, function(x) {x * wcca})
#average for each subject
#sub_mean <- apply(feature_weighted, 1, mean)

#################################################
########### the heatmap without parcellation
#################################################
# methods from Amico & Goni, 2018
##################################

group_connect <- colSums(new_weightMat) #+colMeans(feature_train)
m_connect <- diag(352)
#group_connect <- feature_train[5,]

upperm <- upper.tri(m_connect, diag = F)
m_connect[upperm] <- group_connect
m_connect[lower.tri(m_connect, diag = F)] <- t(m_connect)[lower.tri(m_connect)]
diag(m_connect) <- max(group_connect)

pheatmap(m_connect, treeheight_row = 0, treeheight_col = 0)


#################################################
########### the heatmap with parcellation
#################################################

Nparcels <- sapply(unique(labels_group), function(x) {sum(labels_group == x)})

colnames(m_connect) <- labels_group
rownames(m_connect) <- labels_group
parcels <- unique(labels_group)


# calculate within and between modules connectivity

fig_df <- lapply(seq_along(parcels), function(i){
  
  idxr <- which(rownames(m_connect) == parcels[i])
  M <- sum(rownames(m_connect) == parcels[i])
  within_con <- m_connect[idxr, idxr]
  
  if (i == 21) { # i = 21, brain stem is a special case, only one parcel
    within_load <- within_con
  } else {
    within_load <- 2*sum(within_con)/(M * (M-1))
  }
  
  between_con <- lapply((i+1):length(parcels), function(j) {
    idxc <-  which(colnames(m_connect) == parcels[j])
    between <- m_connect[idxr, idxc]
    M <- sum(rownames(m_connect) == parcels[i])
    N <- sum(colnames(m_connect) == parcels[j])
    list(M = M, N = N, between_con = between, from = parcels[i], to = parcels[j]) 
  })
  
  if (i == length(parcels)) { # cerebellum special case for the last one
    between_load <- NULL
  } else {
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
  }
  df1 <- data.frame(loadings = within_load, from = parcels[i], to = parcels[i])
  df2 <- rbind(df1, between_load)
  list(within_load = within_load, between_load = between_load, fig = df2)
  
})

fig_temp <- do.call(rbind, lapply(fig_df, function(x) {x$fig}))
fig_connecto <- as.data.frame(fig_temp)


within <- do.call(rbind, lapply(fig_df, function(x) {x$within_load}))
between <- do.call(rbind, lapply(fig_df, function(x) {x$between_load}))

######################################
###### connectogram
######################################
library("circlize")
library("tidyverse")

grid.col <- c("#F23F3F","#F4A553","#FFD966","#6AA84F","#42DDC9","#93C6F6","#B39BF7",
              "#EA7BB5","#62E5A0","#F1C232","#F6AC1C","#F6727E","#B0A2DB","#93C6F6",
              "#B7E73C","#F2F618","#F9CB9C","#EA9999","#D5A6BD","#F23F3F","#F4A553",
              "#F1C232","#F6AC1C")



group_fig[group_fig < 1.8 & group_fig > -1.8] <- 0

col_fun <-  colorRamp2(range(group_fig), c("#00BFFF", "#FFD700"), transparency = 0.5)

chordDiagram(group_fig, col=col_fun, grid.col = grid.col)


######################################
###### heatmap
######################################
library('RColorBrewer')
group_connect <- matrix(0,length(parcels),length(parcels))
colnames(group_connect) <- parcels
rownames(group_connect) <- parcels

group_connect[lower.tri(group_connect)] <- -scale(between$loadings)
group_fig <- t(group_connect)
group_fig[lower.tri(group_fig)] <- -scale(between$loadings)
diag(group_fig) <- -scale(within)
#diag(group_fig) <- scale(within, scale = TRUE)
breaksList <- seq(-2, 3, by = 0.01)

pheatmap(group_fig, treeheight_row = 0, treeheight_col = 0, 
         cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)






