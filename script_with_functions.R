### Cleaned data

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# read in files
jornada <- read.csv("data/jornada_rodents.csv", header = TRUE, na.strings = ".")
head(jornada)

sevilleta <- read.csv("data/sevilleta_sm_mrc.txt", header = TRUE)
head(sevilleta)

hj_andrews <- read.csv("data/hj_andrews.csv", header = TRUE)
head(hj_andrews)

shortgrass_steppe <- read.table("data/SGS_LTER_smammals.txt", header = TRUE, sep = "\t", na.strings = ".")
head(shortgrass_steppe)

# prep files

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

# write function for correct grouping of data

group_data