### Cleaned data
# EKB
# March 7, 2016

###################
# LIBRARIES

library(dplyr)
library(tidyr)
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

select_data <- function(data){
  # function for selecting only relevant columns
  dat <- select(data, year, plot, subplot, species)
  return(dat)
}

jor_data <- select_data(jor_data)
hja_data <- select_data(hja_data)
sev_data <- select_data(sev_data)
sgs_data <- select_data(sgs_data)

# list of datasets
datasets <- list(jor_data, sev_data, hja_data, sgs_data)

#####################
# FUNCTIONS

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

find_years_overall <- function(data){
  # create column for relative number of years a species occurs overall
  dat <- select(data, species, year) %>% 
    group_by(species, year) %>% 
    summarise(year_count = n())
  dat1 <- select(dat, species) %>% 
    group_by(species) %>% 
    summarise(years = n())
  total_years <- length(unique(data$year))
  dat1 <- mutate(dat1, rel_years = years/total_years)
  return(dat1)
}

find_locs_overall <- function(data){
  # create column for relative number of locations a species occurs overall
  dat <- select(data, species, siteID) %>% 
    group_by(species, siteID) %>% 
    summarise(site_count = n())
  dat1 <- select(dat, species) %>% 
    group_by(species) %>% 
    summarise(sites = n())
  total_sites <- length(unique(data$siteID))
  dat1 <- mutate(dat1, rel_sites = sites/total_sites)
  return(dat1)
}

join_year_loc <- function(data){
  # combine the outputs of find_years_overall and find_locs_overall into one data frame
  dat <- left_join(x = find_years_overall(data), y = find_locs_overall(data), by = "species")
  return(dat)
}


#####################
# APPLY FUNCTIONS TO DATASETS

jor_data <- group_data(jor_data)
jor_data <- add_siteID(jor_data)
jor_data <- join_year_loc(jor_data)

###################################################################
# WORKING AREA




