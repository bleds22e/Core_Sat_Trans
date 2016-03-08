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

# prep files to make them all similar

jor_data <- select(jornada, year, season, habitat, web, spp, recap) %>% 
            filter(recap != "Y", spp != "DIPO1", spp != "PERO1", spp != "NA") %>% 
            group_by(spp, year, season, habitat, web) %>% 
            arrange(spp, year, season, habitat, web)

