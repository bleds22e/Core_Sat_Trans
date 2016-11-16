# Clustering Algorithms
# Sept 29, 2016
# EKB

# run 'abundances.r'
# use all_abund_data

#####################
# PACKAGES

library(tidyr)
library(dplyr)
library(raster)
library(cluster)
library(vegan)
library(pvclust)

# combine LTER and species columns
all_data <- unite(all_data, sp_LTER, species, LTER, sep = " ")
cluster_data <- all_data %>% dplyr::select(sp_LTER, rel_reg_persist, rel_reg_occup)
