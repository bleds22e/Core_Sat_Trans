# Calculating Abundances
# EKB
# Sept 14, 2016

###################
# LIBRARIES

library(tidyr)
library(dplyr)
library(ggplot2)

###################
# LOAD FILES

jornada <- read.csv("data/Rodents/jornada_rodents.csv", header = TRUE, na.strings = ".")
sevilleta <- read.csv("data/Rodents/sevilleta_sm_mrc.txt", header = TRUE)
hj_andrews <- read.csv("data/Rodents/hj_andrews.csv", header = TRUE)
shortgrass_steppe <- read.table("data/Rodents/SGS_LTER_smammals.txt", header = TRUE, sep = "\t", na.strings = ".")

####################
# PREP FILES

jor_data <- select(jornada, year, season, habitat, web, spp, recap) %>% 
  filter(recap != "Y", spp != "DIPO1", spp != "PERO1", spp != "NA")

hja_data <- select(hj_andrews, YEAR, REPBLK, SUBPLOT, SPECIES) %>% 
  filter(SPECIES != 'UNKN') 

sev_data <- select(sevilleta, year, season, location, web, species, recap) %>% 
  filter(recap != "y", species != "pgsp", species != "dipo", 
         species != "nesp", species != "onsp", species != "pesp", 
         species != "resp", species != "na", species != "pmsp", species != "spsp", 
         year < 2009, location != "two22", location != "blugrama", 
         location != "savanna", location != "goatdraw", location != "rsgrass")

sgs_data <- select(shortgrass_steppe, YEAR, VEG, WEB, SPP) %>% 
  filter(SPP != 'NA')

# rename columns for consistency

names(jor_data) <- c("year", "season", "plot", "subplot", "species", "recap")
names(hja_data) <- c("year", "plot", "subplot", "species")
names(sev_data) <- c("year", "season", "plot", "subplot", "species", "recap")
names(sgs_data) <- c("year", "plot", "subplot", "species")

#####################
# FUNCTIONS

# DATA functions

select_data <- function(data){
  # function for selecting only relevant columns
  dat <- select(data, year, plot, subplot, species)
  return(dat)
}

group_data <- function(data){
  # group and arrange a dataset by species, year, plot, and subplot
  dat <- data %>% 
    group_by(species, year, plot, subplot) %>% 
    arrange(species, year, plot, subplot)
  return(dat)
}

add_siteID <- function(data){
  # make a column with a unique location code
  dat <- data %>% mutate(siteID = paste(plot, subplot, sep = "_"))
  return(dat)
}

prep_data <- function(data){
  # run select_data, group_data, and add_siteID
  data1 <- select_data(data)
  data2 <- group_data(data1)
  data3 <- add_siteID(data2)
  return(data3)
}

adj_abundance <- function(data){
  # calculate the abundance per species, adjusted by the largest average
  abund <- data %>% select(year, siteID, species) %>% 
         group_by(species) %>% 
         summarise(count = n()) %>% 
         mutate(adj_abund = count/max(count))
  return(abund)
}

adj_avg_abund <- function(data){
  # calculate average abundance per species by site/year
  # adjusted by the largest average
  site_year_abund <- data %>% select(year, siteID, species) %>% 
             group_by(species, year, siteID) %>% 
             summarise(count = n())
  avg_abund <- site_year_abund %>% 
             group_by(species) %>% 
             summarise(avg_abund = mean(count)) %>% 
             mutate(ajd_avg_abund = avg_abund/max(avg_abund))
  return(avg_abund)
}

####################
# ABUNDANCES

jor_data <- prep_data(jor_data)
hja_data <- prep_data(hja_data)
sgs_data <- prep_data(sgs_data)
sev_data <- prep_data(sev_data)

# total overall abundance
jor_data_abund <- adj_abundance(jor_data)
hja_data_abund <- adj_abundance(hja_data)
sgs_data_abund <- adj_abundance(sgs_data)
sev_data_abund <- adj_abundance(sev_data)

# average abund by site/year
jor_data_avg_abund <- adj_avg_abund(jor_data)
hja_data_avg_abund <- adj_avg_abund(hja_data)
sgs_data_avg_abund <- adj_avg_abund(sgs_data)
sev_data_avg_abund <- adj_avg_abund(sev_data)