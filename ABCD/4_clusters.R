##############################################################################
############## ABCD: clustering using brain canonical variate scores #########
##############################################################################


###############################################################
############## brain canonical weights across 10 splits 
###############################################################

####### reorder the canonical variates
rs_train_test_abcd <- readRDS("rs_train_test_abcd.rds")
# the reference split: choose a split as a reference split
load_std <- res.abcd.example$v # the loadings for CBCL

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

##### average brain loadings across 10 splits

b <- lapply(seq_along(new_cv_reorder), function(x) {
  loadMat <- rs_train_test_abcd[[x]]$abcd.train$u[, new_cv_reorder[[x]]$idx_reorder]})

brain_mean <- do.call(cbind, lapply(1:6, 
                                function(i) {
                                  rowMeans(abs(do.call(cbind, lapply(b, function(x) {x[, i]})
                                  )))}))
##### calculate the CV scores
# brain_whole is the original brain connecitivity features of the whole ABCD sample
brain_cvscores_average10 <- scale(brain_whole) %*% brain_mean



############################################
######### clustering 
#############################################

scca_brain <- readRDS("brain_cvscores_average10.rds")
cca_brain_1 <- scca_brain[, 1:3]

# 1. hiearachical clustering
d <- dist(cca_brain_1, method ="euclidean")
res.hc <- hclust(d, method = "ward.D")
clusters <- cutree(res.hc, k = 2)

df.clu.figure <- data.frame(CV1 = unlist(cca_brain_1[, 1]), 
                            CV2 = unlist(cca_brain_1[, 2]),
                            CV3 = unlist(cca_brain_1[, 3]))
df.clu.figure$clusters <- clusters
ggplot(df.clu.figure) + geom_point(aes(x = CV1, y = CV2), colour= clusters) +
  theme_bw()

# Dendrogram
d <- as.dendrogram(res.hc)
d <- d %>% color_branches(k=2) %>% color_labels
plot(d, main = "Dendrogram for hierarchical clustering")


# decide how many clusters
cca_rs_data <- cca_brain_1
hcfit_ch <- NbClust(cca_rs_data, method="ward.D",
                    min.nc = 1, max.nc = 6, index = "ch")
hcfit_sl <- NbClust(cca_rs_data, method="ward.D",
                    min.nc = 1, max.nc = 6, index = "silhouette")


fit_ch <- data.frame(number = 1:6, values = hcfit_ch$All.index)
fit_sl <- data.frame(number = 1:6, values = hcfit_sl$All.index)

# figures
ch <- ggplot(fit_ch, aes(number, values))+ 
  geom_line(color = "steelblue2")+ geom_point(color="navyblue" )+
  labs( x= "Number of clusters", y="variance ratio criterion") +
  ggtitle("Variance ratio criterion (Calinski-Harabasz index)") +
  theme_bw()
sl <- ggplot(fit_sl, aes(number, values))+ 
  geom_line(color = "hotpink1")+ geom_point(color="hotpink3" )+
  labs( x= "Number of clusters", y="Silhouette") +
  ggtitle("Silhoutte") +
  theme_bw()

grid.arrange(ch, sl, ncol=2)

# 2. Fuzzy C means clustre analysislibrary(cluster)

cluster_fuzzy_train <- cmeans(cca_brain_1 , centers=2, iter.max = 100, verbose = FALSE,
                              dist = "euclidean", method = "cmeans", m = 2,
                              rate.par = NULL, weights = 1, control = list())
cluster.df <- data.frame(Hclu=clusters, Fclu_1=cluster_fuzzy_train$membership[,1], 
                         Fclu_2 = cluster_fuzzy_train$membership[,2])
dim(cluster.df)


clu_stable <- cluster.df %>% filter(Fclu_1 > 0.7 | Fclu_2 > 0.7)
dim(clu_stable)

idx_clu <- which(cluster.df$Fclu_1 > 0.7 | cluster.df$Fclu_2 > 0.7)

RS_stable <- data.frame(CV1 = unlist(cca_brain_1[-idx,1]), CV2 = unlist(cca_brain_1[-idx,2]),
                        CV3 = unlist(cca_brain_1[-idx,3]))
RS_stable$clusters <- as.factor(clusters[-idx])

idx <- which(cca_brain_1[,2]== min(cca_brain_1[,2]))



plot_ly(RS_stable, x = ~CV1, y = ~CV2, z = ~CV3, 
        color = ~clusters, marker = list(size=3.5), colors = c("skyblue2", "gold1"))

dim(cca_rs_data)
# significance test
sigma <- cov(cca_rs_data)
mu <- colMeans(cca_rs_data)
real_CI <- cluster_test(cca_rs_data)

cl <- makePSOCKcluster(8)
registerDoParallel(cl)
n_sims <- 99
null <- foreach(i = 1:n_sims, .packages=c("NbClust","MASS")) %dopar% {
  rand_sample <- mvrnorm(n=nrow(cca_rs_data), mu=mu, Sigma=sigma)
  res <- cluster_test(rand_sample)
}

stopCluster(cl)



cluster_test <- function(cca_data){
  #ugly hack, because i don't know how to prevent this library creating many plots
  hcfit <- NbClust(cca_data, method="ward.D", index="ch", min.nc=1, max.nc = 6)
  CH_index <- max(hcfit$All.index, na.rm = TRUE)
  hcfit <- NbClust(cca_data, method="ward.D", index="silhouette", min.nc=1, max.nc = 6)
  sil_index <- max(hcfit$All.index, na.rm = TRUE)
  return(c("CH"=CH_index, "Silhouette"=sil_index))
}






null_CI <- as.data.frame(do.call(rbind, null))

rank_cv1 <- sum(real_CI[1] < null_CI[,1]) + 1
pval_cv1 <- rank_cv1 / (n_sims+1)
rank_cv2 <- sum(real_CI[2] < null_CI[,2]) + 1
pval_cv2 <- rank_cv2 / (n_sims+1)
t(t((c("p.val variance ratio"=pval_cv1, "p.val Silhouette"=pval_cv2))))


df.ch <- data.frame(cor = null_CI[,1])
df.sl <- data.frame(cor = null_CI[,2])

ch1 <- ggplot(df.ch, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=20, fill = "gold1", color = "black", alpha = 0.5)+ 
  geom_vline(xintercept=real_CI[1], linetype="dashed", color = "red", size = 1) +
  ggtitle("Variance ratio criterion null") + 
  labs(x = "null distribution") + 
  annotate(geom="text", x=1300, y=60, label="p < 0.001", fontface = 'italic', size = 5, color="black") + 
  theme_bw()
sl2 <- ggplot(df.sl, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.01, fill = "lightblue", color = "black", alpha = 0.5)+
  geom_vline(xintercept=real_CI[2], linetype="dashed", color = "red", size = 1) + 
  annotate(geom="text", x=0.6, y=145, label="p < 0.001", fontface = 'italic', size = 5, color="black") +
  ggtitle("Silhouette null") + 
  labs(x = "null distribution") + 
  theme_bw()

grid.arrange(ch1, sl2, ncol =2)




###############################################
######### clustering: demogrpahic distributions 
###############################################
setwd("~/Desktop/ABCD_download/data/newdata_ReleaseQC")

all_demo <- readRDS("all_final_noNA_incidental_6529.rds")
all_demo$clusters <- clusters

demo_cbcl <- all_demo[,c("sex","race_ethnicity","parental_education",
                        "cbcl_scr_syn_internal_r","cbcl_scr_syn_external_r",
                        "cbcl_scr_syn_totprob_r","clusters")]
summary(demo_cbcl[demo_cbcl$clusters == 1,])
lapply(demo_cbcl[demo_cbcl$clusters == 1,4:6], sd)





