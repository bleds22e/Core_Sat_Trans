# Occupancy vs Time Plots
# Ellen Bledsoe
# April 27, 2016

###################
# LIBRARIES

library(tidyr)
library(dplyr)
library(ggplot2)

###################
# LOAD FILES

jornada <- read.csv("data/jornada_rodents.csv", header = TRUE, na.strings = ".")
sevilleta <- read.csv("data/sevilleta_sm_mrc.txt", header = TRUE)
hj_andrews <- read.csv("data/hj_andrews.csv", header = TRUE)
shortgrass_steppe <- read.table("data/SGS_LTER_smammals.txt", header = TRUE, sep = "\t", na.strings = ".")

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
         location != "savanna", location != "goatdraw")

sgs_data <- select(shortgrass_steppe, YEAR, VEG, WEB, SPP) %>% 
  filter(SPP != 'NA')

# rename columns for consistency

names(jor_data) <- c("year", "season", "plot", "subplot", "species", "recap")
names(hja_data) <- c("year", "plot", "subplot", "species")
names(sev_data) <- c("year", "season", "plot", "subplot", "species", "recap")
names(sgs_data) <- c("year", "plot", "subplot", "species")

####################
# CALCULATIONS

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

occupancy <- function(data){
  # calculate occupancy of species per year, including zeros
  occ <- select(data, year, species, siteID) %>% 
    group_by(species, year) %>% 
    summarize(occupancy = n_distinct(siteID))
  sp <- unique(occ$species)
  yr <- unique(occ$year)
  full_grid <- expand.grid(species = sp, year = yr)
  occ1 <- right_join(occ, full_grid, by = c("species", "year"))
  occ1$occupancy[is.na(occ1$occupancy)] <- 0
  return(occ1)
}


hja_data <- select_data(hja_data)
hja_data <- group_data(hja_data)
hja_data <- add_siteID(hja_data)

#####################
# WORKING



