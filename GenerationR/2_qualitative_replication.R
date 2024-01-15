####################################################################################
################### Greneration R: qualitative repplication ########################
####################################################################################
# This is the total sample of generaiton R
# we ran the whole analysis pipeline here


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
### read the data: brain PCs and CBCL data
brain_genr <- readRDS("brain_genr.rds")
dim(brain_genr)
cbcl_genr <- readRDS("cbcl_genr.rds")
dim(cbcl_genr)

###################################
####### 2. Fit the sCCA model  
###################################
# the penalty parameters were selected based on the most selected parameters in ABCD 

res.genr <- CCA(x=brain_genr, z=cbcl_genr, penaltyx = 0.7, penaltyz = 0.5, typex="standard", typez="standard",
                niter = 20, K=8)

# visualize the loadings:
cbcl_loading <- abs(res.genr$v)
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
                                    nperm=1999, 0.7, 0.5, 8, res.wei.syn$cors)

##################################
######### circular radar plot
##################################
temp <- abs(res.genr$v)
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
######### calculate cross-cohort correlations between ABCD and Generation R CBCL variate scores
##################################

cor_cbcl_scores <- function(x,y, nperm){
    # x: cbcl loadings from GenR
    # y: averaged cbcl loadings from ABCD
    # from ABCD to GenR
    score1 <- as.vector(x %*%t(cbcl_genr_rsfmri))
    score2  <- as.vector(y %*%t(cbcl_genr_rsfmri)) 
    corr_genr <- cor(score1,score2)
    r.per <- replicate(nperm, expr = cor (x = score1, y = sample (score2)))
    p_genr <- sum(abs(r.per) >= abs(corr_genr))/(nperm+1)
    # from GenR to ABCD
    score1 <- as.vector(x %*%t(cbcl_abcd_rsfmri))
    score2 <- as.vector(y %*%t(cbcl_abcd_rsfmri))
    corr_abcd <- cor(score1,score2)
    r.per <- replicate(nperm, expr = cor (x = score1, y = sample (score2)))
    p_abcd <- sum(abs(r.per) >= abs(corr_abcd))/(nperm+1)
    return(list(cor_genr=corr_genr, p_genr=p_genr, cor_abcd=corr_abcd, p_abcd=p_abcd))
}

### attention problems
nperm <- 4999
# from ABCD to GenR
# c_mean is the mean cbcl loadings from ABCD
cor_att <- cor_cbcl_scores(res.genr$v[,4], c_mean[,1], nperm)

### aggressive/rule breaking behaviors
cor_agg <- cor_cbcl_scores(res.genr$v[,2], c_mean[,2], nperm)

### withdrawn behaviors
cor_agg <- cor_cbcl_scores(res.genr$v[,3], c_mean[,3], nperm)

