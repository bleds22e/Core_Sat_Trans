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
  dat <- select(data, year, plot, subplot, species, weight)
  return(dat)
}

avg_weight <- function(data){
  dat <- select(data, species, weight) %>% 
         group_by(species) %>% 
         summarise(avg_weight = mean(weight))
  return(dat)
}

add_systemID <- function(data){
  data$sp_system <- tidyr::unite(species, LTER, sep = "_")
}

add_weight_system <- function(data){
  dat <- left_join(x = data, y = avg_weight(data))
  dat <- add_systemID(data)
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
  dat <- left_join(x = combine_local_reg(data), y = join_abund(data), by = "species")
}

combine_all <- function(data){
  dat <- left_join(x = add_abund(data), y = add_weight_system(data), by = "species")
}

all_together <- function(data){
  # run all functions together to get the full output
  dat <- prep_data(data) %>% combine_all()
  return(dat)
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

systems <- list(c(jor_data, hja_data, sev_data, sgs_data))
full_df <- data.frame()

for (data in systems){
  dat <- all_together(data)
  full_df <- append(dat)
}


####################
# WORK AREA

library(tidyr)
library(dplyr)


# get occupancy and persistence (both local and regional)
jor_data2 <- all_together(jor_data) 
sev_data2 <- all_together(sev_data)
hja_data2 <- all_together(hja_data)
sgs_data2 <- all_together(sgs_data)



# add a dataset column
jor_data2$LTER <- "Jornada Basin"  
hja_data2$LTER <- "H.J. Andrews"
sev_data2$LTER <- "Sevilleta"
sgs_data2$LTER <- "Shortgrass Steppe"

# prep data for calculating abundances
jor_data1 <- prep_data(jor_data)
hja_data1 <- prep_data(hja_data)
sgs_data1 <- prep_data(sgs_data)
sev_data1 <- prep_data(sev_data)
jor_data1$LTER <- "jor"  
hja_data1$LTER <- "hja"
sev_data1$LTER <- "sev"
sgs_data1$LTER <- "sgs"

# total overall abundance
jor_data_abund <- adj_abundance(jor_data1)
hja_data_abund <- adj_abundance(hja_data1)
sgs_data_abund <- adj_abundance(sgs_data1)
sev_data_abund <- adj_abundance(sev_data1)

# average abund by site/year
jor_data_avg_abund <- adj_avg_abund(jor_data1)
hja_data_avg_abund <- adj_avg_abund(hja_data1)
sgs_data_avg_abund <- adj_avg_abund(sgs_data1)
sev_data_avg_abund <- adj_avg_abund(sev_data1)

# add abundance columns to dfs with persistence and occupancy
jor_abund <- left_join(jor_data2, jor_data_abund, by = "species")
jor_abund <- left_join(jor_abund, jor_data_avg_abund, by = "species")

hja_abund <- left_join(hja_data2, hja_data_abund, by = "species")
hja_abund <- left_join(hja_abund, hja_data_avg_abund, by = "species")

sgs_abund <- left_join(sgs_data2, sgs_data_abund, by = "species")
sgs_abund <- left_join(sgs_abund, sgs_data_avg_abund, by = "species")

sev_abund <- left_join(sev_data2, sev_data_abund, by = "species")
sev_abund <- left_join(sev_abund, sev_data_avg_abund, by = "species")

all_abund_data <- bind_rows(jor_abund, hja_abund, sgs_abund, sev_abund)

###################
# WEIGHT and SYSTEM ID



####################
# K-MEANS CLUSTERING

abund_data <- tidyr::unite(all_abund_data, col = sp_system, species, LTER, sep = "_")
abund_data <- tibble::column_to_rownames(all_abund_data, var = "sp_system")
z_data <- scores(all_abund_data)
