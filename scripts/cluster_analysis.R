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
  dat <- ungroup(dat, year, plot, subplot) %>% select(species, sp_system) %>% unique() 
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
  dat <- left_join(x = data, y = add_weight_system(data), by = "species")
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

###################
# GET ABUNDANCES

jor_data <- all_together(jor_data)
hja_data <- all_together(hja_data)
sev_data <- all_together(sev_data)
sgs_data <- all_together(sgs_data)

###################
# ADD WEIGHT AND SYSTEM ID

jor_data$LTER <- "jor"  
hja_data$LTER <- "hja"
sev_data$LTER <- "sev"
sgs_data$LTER <- "sgs"

jor_data <- combine_all(jor_data)
hja_data <- combine_all(hja_data)
sev_data <- combine_all(sev_data)
sgs_data <- combine_all(sgs_data)

####################
# WORK AREA


# test <- all_together(jor_data)
test <- prep_data(jor_data)
test <- combine_all(jor_data)



systems <- list(c(jor_data, hja_data, sev_data, sgs_data))
full_df <- data.frame()

for (i in 1:length(systems)){
  dat <- all_together(i)
  full_df[i,] <- dat
  return(dat)
}

####################
# K-MEANS CLUSTERING

abund_data <- tibble::column_to_rownames(all_abund_data, var = "sp_system")
z_data <- scores(all_abund_data)
