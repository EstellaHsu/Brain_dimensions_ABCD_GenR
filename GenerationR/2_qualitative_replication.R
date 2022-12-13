####################################################################################
################### Greneration R: qualitative repplication ########################
####################################################################################
# This is the total sample of generaiton R
# we will run the whole analysis pipeline here


########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'fmsb', 'NbClust','stats','dendextend','cluster',
              'fpc', 'e1071', 'plotly', 'MASS', 'ggradar2','tidyr','introdatavizyes','Hmisc',
              'pheatmap','RColorBrewer','circlize')


lapply(packages, require, character.only = TRUE)


###################################
####### 1. Read the data 
###################################
### read the data: because I did different imputation based on syndrome scale and item level, I have different residual matrix
brain_genr <- readRDS("brain_genr.rds")
dim(brain_genr)
cbcl_genr <- readRDS("cbcl_genr.rds")
dim(cbcl_genr)

###################################
####### 2. Fit the sCCA model  
###################################
# the penalty parameters were selected based on the most selected parameters in ABCD (see Table 1)

res.genr <- CCA(x=brain_genr, z=cbcl_genr, penaltyx = 0.5, penaltyz = 0.5, typex="standard", typez="standard",
                niter = 20, K=8)

# visualize the loadings:
cbcl_loading <- abs(res.wei.syn$v)
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")

corrplot(t(cbcl_loading)[1:8,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 2, number.cex=1)


###################################
######## Covariance explained  
###################################
vardf <- VarianceExplain(brain_genr, cbcl_genr, res.wei.syn, 8) 
colnames(vardf) <- c("principal_components", "Covariance_explained")
ggplot(vardf, aes(x=principal_components, y=Covariance_explained))+
  geom_point(size = 3) + theme_bw()

###################################
######## permutation test  
###################################

perm_genr_total <- permutation_test(cbcl_genr, brain_genr, 
                                    nperm=1999, 0.5, 0.5, 8, res.wei.syn$cors)

##################################
######### circular radar plot
##################################
temp <- abs(res.wei.syn$v)
df_temp1 <- as.data.frame(temp)
df_temp1 <- t(df_temp1)
df_temp2 <- rbind(rep(1,8), df_temp1)
colnames(df_temp2) <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule breaking","Aggression")

df_temp3 <- as.data.frame(cbind(group=c("1","CV1", "CV2", "CV3", "CV4", "CV5"),df_temp2))
df_temp3[,2:9] <- lapply(df_temp3[,2:9], function(x){as.numeric(as.character(x))})

colors_border_6 <- c("#fcfcfb","#CC99C9","#F7EA48","#4FC1E8","#F08080","#A0D568")
colors_in_6 <- c("#fcfcfb","#CC99C9","#F7EA48","#4FC1E8","#F08080","#A0D568")

p.genr <- ggradar2(df_temp3, base.size=5, webtype = 'lux', grid.min = 0,
               grid.max = 1,label.gridline.mid = TRUE, group.colours = colors_border_6,
               group.fill.colours = colors_in_6,label.centre.y=FALSE,
               gridline.mid.colour="grey", grid.label.size = 0,
               gridline.max.linetype = "solid",polygonfill.transparency=0.5,
               background.circle.transparency=0.1,axis.label.size=5.5,
               group.line.width=1.5,group.point.size=3,plot.legend = TRUE) 
p.genr

ggsave("radarfigure_genr.pdf", width = 12, height = 6)
dev.off()

##################################
######### calculate correlations between ABCD and Generation R CBCL loadings
##################################
for(i in 1:6){
  for(j in 1:5){
    cor_abcd_genr <- cor(cbcl_mean[,i], res.wei.syn$v[,j]) # cbcl_mean is the average CBCL loadings across 10 train-test splits
    out <- paste0("abcd_CV",i,"  genr_CV",j,":", round(cor_abcd_genr,2))
    print(out)
  }
}

