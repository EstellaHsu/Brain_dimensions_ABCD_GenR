########################################################
############## Greneration R ########################
########################################################
# This is the total sample of generaiton R
# we will run the whole analysis pipeline here
# the correct weighted PCA will be used here. 

########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'stats', 'cluster','fpc','e1071','dendextend','NbClust',
              'plotly')

lapply(packages, require, character.only = TRUE)


########################################################
############## 1. Read the data ########################
########################################################
setwd("/Users/estella/Desktop/ABCD_download/data/Generation_R")
### read the data: because I did different imputation based on syndrome scale and item level, I have different residual matrix
feature_genr_syn <- readRDS("feature_residual_aromaonly_cbclsyndrome.rds")
dim(feature_genr_syn)
cbcl_syndrome_genr <- readRDS("GenR_confounders_imp_cbclsyndrome.rds")
names(cbcl_syndrome_genr)
summary(cbcl_syndrome_genr)
table(cbcl_syndrome_genr$gender)
sd(cbcl_syndrome_genr$age)
# select syndrome scales
cbcl_syn_genr <- cbcl_syndrome_genr[,1:8]
names(cbcl_syn_genr)
str(cbcl_syn_genr)
########################################################
############## 2. Weighted PCA  ########################
########################################################
# This is based on the method Hao gave me. 
rank <- rank(-rowSums(cbcl_syn_genr))
weights <- log(nrow(cbcl_syn_genr)) - log(rank)
weights <- weights / sum(weights)
# center the original data
feature_brain_centered <- scale(feature_genr_syn, scale = FALSE)
# scale the columns of the feature matrix
feature_centered_weighted <- feature_brain_centered * replicate(ncol(feature_brain_centered), weights)
# call PCA on the weighted feature matrix
pca.weighted <- prcomp(feature_centered_weighted)
# the rotation transformation: take the first 100 eigenvectors
rotation <- pca.weighted$rotation[, 1:100]
# the data with the reduced dimensionality
feature_brain_reduced <- feature_brain_centered %*% rotation 
f <- feature_brain_centered %*% rotation 
dim(feature_brain_reduced)
vars <- apply(f, 2, var)
cumsum(vars/sum(vars))
#pca.norm <- prcomp(feature_genr_syn)
#var_explained <- pca.norm$sdev^2/sum(pca.norm$sdev^2)
#cumsum(var_explained)[1:200]
pca.weighted <- readRDS("pca.weighted.genr.rds")
saveRDS(rotation, "rotation.rds")
saveRDS(pca.weighted, "pca.weighted.genr.rds")
saveRDS(pca.norm, "pca.norm.genr.rds")
saveRDS(feature_brain_centered, "feature_brain_centered_genr.rds")
pca.weighted <- readRDS("pca.weighted.genr.rds")
feature_brain_centered <- readRDS("feature_brain_centered_genr.rds")

########## or sparce PCA: too slow 
#cv.out <- SPC.cv(feature_genr_syn, sumabsvs=seq(1.2, sqrt(ncol(feature_genr_syn)), len=10), nfolds = 5)
#print(cv.out)
#plot(cv.out)
#out <- SPC(x,sumabsv=cv.out$bestsumabs, K=4) # could use

feature_genr_projectFromabcd <- feature_brain_centered %*% rs_rotation[[7]]


########################################################
############## 3. Grid search for penalty parameters  ##
########################################################

####### Step 1: train penalty parameters 
##### weighted PCA
cv.resample <- CV_sampling(cbcl_syn_genr, feature_genr_projectFromabcd, 100, 0.8)
x_pen <- seq(0.1, 1.0, length.out = 10) # brain
y_pen <- seq(0.1, 1.0, length.out = 10)


grid.syn.wei.genr <- grid.search.cor.Testset(cv.resample$brain_train, cv.resample$cbcl_train,
                                              cv.resample$brain_test, cv.resample$cbcl_test,
                                              x_pen, y_pen,nsample=100)

saveRDS(grid.Testset.syn.wei.genr, "grid.cor.Testset.syn.wei.wholesample.rds")
grid.genr <- readRDS("grid.cor.Testset.syn.wei.wholesample.rds")

########################################################
############## 3. Fit the sCCA model  ##################
########################################################
CCA.permute(x=feature_genr_projectFromabcd, z=cbcl_syn_genr,typex="standard", typez="standard")

res.wei.syn <- CCA(x=feature_brain_reduced, z=cbcl_syn_genr, penaltyx = 0.5, penaltyz = 0.5, typex="standard", typez="standard",
                   niter = 20, K=8)

# visualize the loadings:
cbcl_loading <- abs(res.wei.syn$v)
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")

corrplot(t(cbcl_loading)[1:8,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 2, number.cex=1)

# similar, just with flipped CV
#res.norm.syn <- CCA(x=pca.norm$x[,1:100], z=cbcl_syn_genr, penaltyx = 0.3, penaltyz = 0.8, typex="standard", typez="standard",
                   #niter = 20, K=8)

 
########################################################
############## Covariance explained  ###################
########################################################
vardf <- VarianceExplain(feature_brain_reduced, cbcl_syn_genr, res.wei.syn, 8) 
colnames(vardf) <- c("principal_components", "Covariance_explained")
ggplot(vardf, aes(x=principal_components, y=Covariance_explained))+
  geom_point(size = 3) + theme_bw()

########################################################
############## permutation test  #######################
########################################################
perm_genr_total <- permutation_test(cbcl_syn_genr, feature_brain_reduced, 
                                    nperm=1999, 0.5, 0.5, 8, res.wei.syn$cors)
# lot the correlations
df.cor.2 <- data.frame(cor = perm_genr_total$cor.perm[, 2])
df.cor.1 <- data.frame(cor = perm_genr_total$cor.perm[, 1])


cor1 <- ggplot(df.cor.1, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "gold1", color = "black", alpha = 0.6)+
  geom_vline(xintercept=0.31, linetype="dashed", color = "red", size = 1) +
  ggtitle("Permutation correlations of the 1st canonical variate") + 
  labs(x = "canonical correlations") + 
  annotate(geom="text", x=0.3, y=200, label="r = 0.31", size = 5, color="black") + xlim(c(0.1, 0.4))+
  annotate(geom="text", x=0.3, y=180, label="p(fdr) < 0.001", fontface = 'italic', size = 5, color="black") + 
  theme_bw()
cor2 <- ggplot(df.cor.2, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "lightblue", color = "black", alpha = 0.8)+
  geom_vline(xintercept=0.26, linetype="dashed", color = "red", size = 1) + 
  annotate(geom="text", x=0.4, y=250, label="r = 0.26", size = 5, color="black") + xlim(c(0.1, 0.4))+
  annotate(geom="text", x=0.4, y=225, label="p(fdr) < 0.001", fontface = 'italic', size = 5, color="black") +
  ggtitle("Permutation correlations of the 2nd canonical variate") + 
  labs(x = "canonical correlations") + 
  theme_bw()

grid.arrange(cor1, cor2, ncol =2)  

##################################
######### circular radar plot
##################################
temp <- abs(res.wei.syn$v[,c(5,2,3,1,4)])
df_temp1 <- as.data.frame(temp)
df_temp1 <- t(df_temp1)
df_temp2 <- rbind(rep(1,8), df_temp1)
colnames(df_temp2) <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule breaking","Aggression")

df_temp3 <- as.data.frame(cbind(group=c("1","CV5", "CV2", "CV3", "CV1", "CV4"),df_temp2))
df_temp3[,2:9] <- lapply(df_temp3[,2:9], function(x){as.numeric(as.character(x))})


#pdf("radarfigure_train_abcd.pdf", width=12, height=8)
colors_border_6 <- c("#fcfcfb","#CC99C9","#F7EA48","#4FC1E8","#F08080","#A0D568")
colors_in_6 <- c("#fcfcfb","#CC99C9","#F7EA48","#4FC1E8","#F08080","#A0D568")

p2 <- ggradar2(df_temp3, base.size=5, webtype = 'lux', grid.min = 0,
               grid.max = 1,label.gridline.mid = TRUE, group.colours = colors_border_6,
               group.fill.colours = colors_in_6,label.centre.y=FALSE,
               gridline.mid.colour="grey", grid.label.size = 0,
               gridline.max.linetype = "solid",polygonfill.transparency=0.5,
               background.circle.transparency=0.1,axis.label.size=5.5,
               group.line.width=1.5,group.point.size=3,plot.legend = TRUE) 
p2

ggsave("radarfigure_genr.png", width = 12, height = 6)
dev.off()

##################################
######### calculate correlations
##################################
for(i in 1:6){
  for(j in 1:5){
    cor_abcd_genr <- cor(c_mean[,i], temp[,j])
    out <- paste0("abcd_CV",i,"  genr_CV",j,":", round(cor_abcd_genr,2))
    print(out)
  }
}

# abcd CV1, CV5, 0.79
# abcd CV2, CV3, 0.36
# abcd CV3, CV3, 0.67
# abcd CV6, CV1, 0.95


########################################################
############## Replication in ABCD  ####################
########################################################
feature_total_syn <- readRDS("feature_genrMapToabcd.rds")
dim(feature_total_syn)
#saveRDS(test.pca.genr.total, "test.pca.genr.total.rds")

# direct mapping
cor.genrToABCD <- test_project_weights(feature_genr_projectFromabcd, cbcl_syn_genr, res.wei.syn.abcd, 6)
cor.genrToABCD
perm_abcdTogenr_total <- permutation_test_testset(cbcl_syn_abcd,feature_total_syn[,1:20],nperm=999, res.wei.syn,cor.genrToABCD)
## Step 5.1: predict the PCs in Generation R





###############################################################
############## Visualization of CBCL syndrome loadings  #######
###############################################################


df_genr <- abs(res.wei.syn$v[, 1:5])
df_genr <- as.data.frame(df_genr)
cbclnames <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")

df_genr <- t(df_genr)
colnames(df_genr) <- cbclnames
rownames(df_genr) <- c("CV1", "CV2", "CV3", "CV4", "CV5")
df_cbcl_genr <- rbind(rep(1,8), rep(0,8), df_genr)
df_cbcl_genr <- as.data.frame(df_cbcl_genr)

colors_border <- c( rgb(0.2,0.7,0.6,0.9), rgb(1,0.8,0.1,0.9), rgb(0.5,0.2,0.9,0.9), rgb(0.9,0.5,0.5,0.9), rgb(0.3,0.5,0.8,0.9))
colors_in <- c( rgb(0.3,0.7,0.6,0.4), rgb(1,0.9,0,0.4), rgb(0.7,0.5,0.9,0.4), rgb(1,0.7,0.6,0.4), rgb(0.6,0.8,1,0.4))

# plot with default options:
radarchart( df_cbcl_genr , axistype=1, 
            #custom polygon
            pcol= colors_border , pfcol=colors_in , plwd=5, plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0.2,1.0,0.2), cglwd=0.8,
            #custom labels
            vlcex=1.2
)

# Add a legend
legend(x=1.5, y=1, legend = rownames(df_cbcl1[-c(1,2),]), 
       bty = "n", pch=20, col=colors_border , text.col = "black", cex=1, pt.cex=2)






