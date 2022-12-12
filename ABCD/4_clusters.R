##############################################################################
############## ABCD: clustering using brain canonical variate scores #########
##############################################################################


###############################################################
############## brain canonical weights across 10 splits 
###############################################################

####### reorder the canonical variates
rs_train_test_abcd <- readRDS("rs_train_test_abcd.rds")
# the reference split: choose a split as a reference split
load_ref <- res.abcd.example$v # the loadings for CBCL

new_cv_reorder  <- lapply(seq_along(rs_train_test_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:6) {
    # calculate the maximum correlation of each CV and reorder, we only care about the first six correlations (or any number)
    cor_cv <- sapply(1:6, function(z) {cor(load_ref[, i], rs_train_test_abcd[[j]]$abcd.train$v[, z])})
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

brain <- lapply(seq_along(new_cv_reorder), function(x) {
  loadMat <- rs_train_test_abcd[[x]]$abcd.train$u[, new_cv_reorder[[x]]$idx_reorder]})

brain_mean <- do.call(cbind, lapply(1:6, 
                                function(i) {
                                  rowMeans(abs(do.call(cbind, lapply(brain, function(x) {x[, i]})
                                  )))}))
##### calculate the CV scores
# brain_whole is the original brain connecitivity features of the whole ABCD sample
# i.e., it is the PCs calculated by weighted PCA based on the whole ABCD sample
brain_cvscores_average10 <- scale(brain_whole) %*% brain_mean



############################################
######### clustering 
#############################################

# choose the first 3 brain CVs
scca_brain <- brain_cvscores_average10[, 1:3]

######### the following code are from Dinga et al., 2019

# hiearachical clustering
d <- dist(scca_brain, method ="euclidean")
res.hc <- hclust(d, method = "ward.D")
clusters <- cutree(res.hc, k = 2)

df.clu.figure <- data.frame(CV1 = unlist(cca_brain_1[, 1]), 
                            CV2 = unlist(cca_brain_1[, 2]),
                            CV3 = unlist(cca_brain_1[, 3]))
df.clu.figure$clusters <- clusters


# Dendrogram
d <- as.dendrogram(res.hc)
d <- d %>% color_branches(k=2) %>% color_labels
plot(d, main = "Dendrogram for hierarchical clustering")


# decide how many clusters
# variance ratio criterion
hcfit_ch <- NbClust(scca_brain, method="ward.D",
                    min.nc = 1, max.nc = 6, index = "ch")
# sihoutte
hcfit_sl <- NbClust(scca_brain, method="ward.D",
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

# 3-D visualization 
df.clu.figure <- data.frame(CV1 = unlist(cca_brain_1[, 1]), 
                            CV2 = unlist(cca_brain_1[, 2]),
                            CV3 = unlist(cca_brain_1[, 3]))
df.clu.figure$clusters <- clusters

plot_ly(df.clu.figure, x = ~CV1, y = ~CV2, z = ~CV3, 
        color = ~clusters, marker = list(size=3.5), colors = c("skyblue2", "gold1"))


# significance test
sigma <- cov(scca_brain)
mu <- colMeans(scca_brain)
real_CI <- cluster_test(scca_brain)

cl <- makePSOCKcluster(8)
registerDoParallel(cl)
n_sims <- 999
null <- foreach(i = 1:n_sims, .packages=c("NbClust","MASS")) %dopar% {
      cluster_test <- function(cca_data){
      hcfit <- NbClust(cca_data, method="ward.D", index="ch", min.nc=1, max.nc = 6)
      CH_index <- max(hcfit$All.index, na.rm = TRUE)
      hcfit <- NbClust(cca_data, method="ward.D", index="silhouette", min.nc=1, max.nc = 6)
      sil_index <- max(hcfit$All.index, na.rm = TRUE)
      return(c("CH"=CH_index, "Silhouette"=sil_index))
     }
  rand_sample <- mvrnorm(n=nrow(cca_rs_data), mu=mu, Sigma=sigma)
  res <- cluster_test(rand_sample)
}

stopCluster(cl)


null_CI <- as.data.frame(do.call(rbind, null))

rank_cv1 <- sum(real_CI[1] < null_CI[,1]) + 1
pval_cv1 <- rank_cv1 / (n_sims+1)
rank_cv2 <- sum(real_CI[2] < null_CI[,2]) + 1
pval_cv2 <- rank_cv2 / (n_sims+1)
t(t((c("p.val variance ratio"=pval_cv1, "p.val Silhouette"=pval_cv2))))

# visualization of the significance test
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


