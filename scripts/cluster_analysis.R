# Cluster Analysis Using Abundance Data
# EKB
# November 2016

library(tidyr)
library(dplyr)
library(ggplot2)

###################
# LOAD ABUNDANCE FUNCTIONS

source("scripts/abundance_functions.R")

# Other functions

select_data <- function(data){
  # function for selecting only relevant columns
  dat <- select(data, year, plot, subplot, species, weight, LTER)
  return(dat)
}

avg_weight <- function(data){
  dat <- data %>% group_by(species) %>% 
         summarise(avg_weight = mean(weight))
  return(dat)
}

add_systemID <- function(data){
  dat <- tidyr::unite(data, sp_system, species, LTER, sep = "_", remove = FALSE)
  dat <- ungroup(dat, year, plot, subplot) %>% select(species, sp_system, LTER) %>% unique() 
  return(dat)
}

add_weight_system <- function(data){
  dat <- left_join(x = add_systemID(data), y = avg_weight(data), by = "species")
  return(dat)
}

join_abund <- function(data){
  # combine the relative regional persistance and regional occupancy into one data frame
  dat <- left_join(x = adj_abundance(data), y = adj_avg_abund(data), by = "species")
  return(dat)
}

join_local_reg <- function(data){
  # combine regional and local persistance and occupancy into one data frame
  dat <- left_join(x = join_regionals(data), y = join_locals(data), by = "species")
  return(dat)
}

add_abund <- function(data){
  dat <- left_join(x = join_local_reg(data), y = join_abund(data), by = "species")
}

all_together <- function(data){
  # run all functions together to get the full output
  dat <- prep_data(data) %>% add_abund()
  return(dat)
}

combine_all <- function(data){
  original <- prep_data(data)
  new <- all_together(data)
  dat <- left_join(x = new, y = add_weight_system(original), by = c("species"))
}

###################
# LOAD FILES

jornada <- read.csv("data/Rodents/jornada_rodents.csv", header = TRUE, na.strings = ".")
sevilleta <- read.csv("data/Rodents/sevilleta_sm_mrc.txt", header = TRUE)
hj_andrews <- read.csv("data/Rodents/hj_andrews.csv", header = TRUE)
shortgrass_steppe <- read.table("data/Rodents/SGS_LTER_smammals.txt", header = TRUE, sep = "\t", na.strings = ".")

####################
# PREP FILES

jor_data <- select(jornada, year, habitat, web, spp, recap, weight) %>% 
  filter(recap != "Y", spp != "DIPO1", spp != "PERO1", spp != "NA")

hja_data <- select(hj_andrews, YEAR, REPBLK, SUBPLOT, SPECIES, WEIGHT) %>% 
  filter(SPECIES != 'UNKN') 

sev_data <- select(sevilleta, year, location, web, species, recap, mass) %>% 
  filter(recap != "y", species != "pgsp", species != "dipo", 
         species != "nesp", species != "onsp", species != "pesp", 
         species != "resp", species != "na", species != "pmsp", species != "spsp", 
         year < 2009, location != "two22", location != "blugrama", 
         location != "savanna", location != "goatdraw", location != "rsgrass")

sgs_data <- select(shortgrass_steppe, YEAR, VEG, WEB, SPP, WT) %>% 
  filter(SPP != 'NA')

# rename columns for consistency

names(jor_data) <- c("year", "plot", "subplot", "species", "recap", "weight")
names(hja_data) <- c("year", "plot", "subplot", "species", "weight")
names(sev_data) <- c("year", "plot", "subplot", "species", "recap", "weight")
names(sgs_data) <- c("year", "plot", "subplot", "species", "weight")

jor_data$LTER <- "jor"  
hja_data$LTER <- "hja"
sev_data$LTER <- "sev"
sgs_data$LTER <- "sgs"

###################
# GET ALL DATA

jor_data <- combine_all(jor_data)
hja_data <- combine_all(hja_data)
sev_data <- combine_all(sev_data)
sgs_data <- combine_all(sgs_data)

# merge into one dataframe

all_data <- bind_rows(jor_data, hja_data, sev_data, sgs_data)

####################
# PREP FOR CLUSTERING

library(raster)
library(cluster)

# add row names

abund_data <- tibble::column_to_rownames(all_data, var = "sp_system") 
abund_data <- dplyr::select(abund_data, rel_reg_persist, rel_reg_occup, mean_local_persist,
                            mean_local_occup, adj_abund, adj_avg_abund)

# check normality

par(mfrow = c(3,2))
hist(all_data$rel_reg_persist)
hist(all_data$rel_reg_occup)
hist(all_data$mean_local_persist)
hist(all_data$mean_local_occup)
hist(all_data$adj_abund)
hist(all_data$adj_avg_abund)

# transform data

#trans1_rel_reg_persist <- log(asin(sqrt(reflect+1)))
#trans1_rel_reg_occup <- log(asin(sqrt(reflect+1)))
#trans_mean_local_persist <- asin(sqrt(all_data$mean_local_persist))
#trans_mean_local_occup <- asin(sqrt(all_data$mean_local_occup))
#trans_adj_abund <- log(all_data$adj_abund)
#trans_adj_avg_abund <- log(all_data$adj_avg_abund)

# standardize data

z_data <- scale(abund_data)
z_regional <- dplyr::select(as.data.frame(z_data), rel_reg_persist, rel_reg_occup, adj_abund)
z_local <- dplyr::select(as.data.frame(z_data), mean_local_persist, mean_local_occup, adj_avg_abund)

## Plots for Number of Clusters
par(mfrow = c(2,2))

# Regional Data

# make a scree plot
scree_regional <- rep(0, 20)
for (i in 1:20){
  scree_regional[i] <- sum(kmeans(z_regional, center = i, nstart = 25)$withinss)
}
plot(1:10, scree_regional[1:10], type = "b", xlab = "Number of groups", 
     ylab = "Within groups sum of squares") 

# silhouette plot
sil_regional <- rep(0,20)
for (i in 2:20){
  sil_regional[i] <- summary(silhouette(kmeans(z_regional, centers = i, iter.max = 100, 
                                      nstart = 25)$cluster, dist(z_data)))$avg.width
}
plot(2:10, sil_regional[2:10], type = "b", xlab = "Number of groups", ylab = "average silhouette width ")

# Local Data

# make a scree plot
scree_local <- rep(0, 20)
for (i in 1:20){
  scree_local[i] <- sum(kmeans(z_local, center = i, nstart = 25)$withinss)
}
plot(1:10, scree_local[1:10], type = "b", xlab = "Number of groups", 
     ylab = "Within groups sum of squares") 

# silhouette plot
sil_local <- rep(0,20)
for (i in 2:20){
  sil_local[i] <- summary(silhouette(kmeans(z_local, centers = i, iter.max = 100, 
                                     nstart = 25)$cluster, dist(z_data)))$avg.width
}
plot(2:10, sil_local[2:10], type = "b", xlab = "Number of groups", ylab = "average silhouette width ")

#######################
# K-MEANS CLUSTERING

### Regional -- 3 CLUSTERS
k_regional <- kmeans(z_regional, centers = 3, iter.max = 10, nstart = 25)
pairs(z_regional, panel = function(x, y, z) text(x, y, k_regional$cluster))

# plot against PC 1&2
  # run PCA
pca_regional <- princomp(z_regional, cor=T)
summary(pca_regional)
pca_regional$loadings
  # colors
my.color.vector <- rep("red", times=nrow(z_regional))
my.color.vector[k_regional$cluster==1] <- "red"
my.color.vector[k_regional$cluster==2] <- "green"
my.color.vector[k_regional$cluster==3] <- "blue"
  # plot clusters
plot(pca_regional$scores[,1], pca_regional$scores[,2], ylim=range(pca_regional$scores[,1]),xlim=range(pca_regional$scores[,1]*1.25), xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
text(pca_regional$scores[,1], pca_regional$scores[,2], labels=rownames(z_regional), cex=1.25, lwd=2,
     col=my.color.vector)
  # plot onto biplot
biplot(pca_regional)
text(pca_regional$scores[,1], pca_regional$scores[,2], labels=rownames(z_regional), cex=1.25, lwd=2,
     col=my.color.vector)

# different way to plot

scores <- pca_regional$scores
  # gg: data frame of PC1 and PC2 scores with corresponding cluster
gg <- data.frame(cluster=factor(k_regional$cluster), x=scores[,1], y=scores[,2])
  # calculate cluster centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
  # merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
  # calculate 95% confidence ellipses
library(ellipse)
conf.rgn  <- do.call(rbind,lapply(1:3,function(i)
  cbind(cluster=i,ellipse(cov(gg[gg$cluster==i,2:3]),centre=as.matrix(centroids[i,2:3])))))
conf.rgn  <- data.frame(conf.rgn)
conf.rgn$cluster <- factor(conf.rgn$cluster)
  # plot cluster map
reg_clust_3_plot <- ggplot(gg, aes(x,y, color=cluster))+
  geom_point(size=3) +
  geom_point(data=centroids, size=4) +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y))+
  geom_path(data=conf.rgn)
reg_clust_3_plot
  # add cluster id to dataframe
reg_clust_3 <- as.data.frame(k_regional$cluster)
names(reg_clust_3) <- c("reg_clust_3")
all_clust_data <- dplyr::bind_cols(all_data, reg_clust_3)

# where the clusters fall on occupancy x persistence plot
ggplot(all_clust_data, aes(x = rel_reg_occup, y = rel_reg_persist)) +
  geom_point(aes(color = as.character(reg_clust_3)), size = 3)

#============================================

### Regional -- 2 CLUSTERS
k2_regional <- kmeans(z_regional, centers = 2, iter.max = 10, nstart = 25)
pairs(z_regional, panel = function(x, y, z) text(x, y, k2_regional$cluster))

# plot against PC 1&2
  # colors
my.color.vector <- rep("red", times=nrow(z_regional))
my.color.vector[k2_regional$cluster==1] <- "red"
my.color.vector[k2_regional$cluster==2] <- "blue"
  # plot clusters
plot(pca_regional$scores[,1], pca_regional$scores[,2], ylim=range(pca_regional$scores[,1]),xlim=range(pca_regional$scores[,1]*1.25), xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
text(pca_regional$scores[,1], pca_regional$scores[,2], labels=rownames(z_regional), cex=1.25, lwd=2,
     col=my.color.vector)
  # plot onto biplot
biplot(pca_regional)
text(pca_regional$scores[,1], pca_regional$scores[,2], labels=rownames(z_regional), cex=1.25, lwd=2,
     col=my.color.vector)

# different way to plot

gg2 <- data.frame(cluster=factor(k2_regional$cluster), x=scores[,1], y=scores[,2])
  # calculate cluster centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg2,mean)
  # merge centroid locations into ggplot dataframe
gg2 <- merge(gg2,centroids,by="cluster",suffixes=c("",".centroid"))
  # calculate 95% confidence ellipses -- NOT WORKING YET
#conf.rgn  <- do.call(rbind,lapply(1:2,function(i)
#  cbind(cluster=i,ellipse(cov(gg2[gg2$cluster==i,2]),centre=as.matrix(centroids[i,2])))))
#conf.rgn  <- data.frame(conf.rgn)
#conf.rgn$cluster <- factor(conf.rgn$cluster)
  # plot cluster map
reg_clust_2_plot <- ggplot(gg2, aes(x,y, color=cluster))+
  geom_point(size=3) +
  geom_point(data=centroids, size=4) +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y))

reg_clust_2_plot
  # add cluster id to dataframe
reg_clust_2 <- as.data.frame(k2_regional$cluster)
names(reg_clust_2) <- c("reg_clust_2")
all_clust_data <- dplyr::bind_cols(all__clust_data, reg_clust_2)

# where the clusters fall on occupancy x persistence plot
ggplot(all_clust_data, aes(x = rel_reg_occup, y = rel_reg_persist)) +
  geom_point(aes(color = as.character(reg_clust_2)), size = 3)

#############################
# AVERAGE MASS BY CLUSTER

reg_by_cluster_3 <- all_clust_data %>% # weight is showing up as NA?
                    dplyr::select(reg_clust_3, rel_reg_persist, rel_reg_occup, adj_abund, avg_weight) %>% 
                    dplyr::group_by(reg_clust_3) %>% 
                    dplyr::summarise_each(funs(mean))
head(reg_by_cluster_3)

#=====================================
z_regional <- dplyr::select(as.data.frame(z_data), rel_reg_persist, rel_reg_occup, adj_avg_abund)

### Regional -- 3 CLUSTERS
k_regional <- kmeans(z_regional, centers = 3, iter.max = 10, nstart = 25)
pairs(z_regional, panel = function(x, y, z) text(x, y, k_regional$cluster))

# plot against PC 1&2
# run PCA
pca_regional <- princomp(z_regional, cor=T)
summary(pca_regional)
pca_regional$loadings
# colors
my.color.vector <- rep("red", times=nrow(z_regional))
my.color.vector[k_regional$cluster==1] <- "red"
my.color.vector[k_regional$cluster==2] <- "green"
my.color.vector[k_regional$cluster==3] <- "blue"
# plot clusters
plot(pca_regional$scores[,1], pca_regional$scores[,2], ylim=range(pca_regional$scores[,1]),xlim=range(pca_regional$scores[,1]*1.25), xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
text(pca_regional$scores[,1], pca_regional$scores[,2], labels=rownames(z_regional), cex=1.25, lwd=2,
     col=my.color.vector)
# plot onto biplot
biplot(pca_regional)
text(pca_regional$scores[,1], pca_regional$scores[,2], labels=rownames(z_regional), cex=1.25, lwd=2,
     col=my.color.vector)

# different way to plot

scores <- pca_regional$scores
# gg: data frame of PC1 and PC2 scores with corresponding cluster
gg <- data.frame(cluster=factor(k_regional$cluster), x=scores[,1], y=scores[,2])
# calculate cluster centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)
# merge centroid locations into ggplot dataframe
gg <- merge(gg,centroids,by="cluster",suffixes=c("",".centroid"))
# calculate 95% confidence ellipses
library(ellipse)
conf.rgn  <- do.call(rbind,lapply(1:3,function(i)
  cbind(cluster=i,ellipse(cov(gg[gg$cluster==i,2:3]),centre=as.matrix(centroids[i,2:3])))))
conf.rgn  <- data.frame(conf.rgn)
conf.rgn$cluster <- factor(conf.rgn$cluster)
# plot cluster map
reg_clust_3_plot <- ggplot(gg, aes(x,y, color=cluster))+
  geom_point(size=3) +
  geom_point(data=centroids, size=4) +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y))+
  geom_path(data=conf.rgn)
reg_clust_3_plot
# add cluster id to dataframe
reg_clust_3 <- as.data.frame(k_regional$cluster)
names(reg_clust_3) <- c("reg_clust_3")
all_clust_data <- dplyr::bind_cols(all_data, reg_clust_3)

# where the clusters fall on occupancy x persistence plot
ggplot(all_clust_data, aes(x = rel_reg_occup, y = rel_reg_persist)) +
  geom_point(aes(color = as.character(reg_clust_3)), size = 3)

