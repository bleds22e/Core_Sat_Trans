# Core-satellite-transient with Portal data
# Feb 17, 2016
# Ellen Bledsoe

# download libraries
library(dplyr)
library(ggplot2)

# read files

species <- read.csv("data/species_portal.csv", header = TRUE, na.strings = c(""))
surveys <- read.csv("data/surveys.csv", header = TRUE, na.strings = c(""))

# get only known rodents and appropriate periods

rodents <- select(species, New.Code, Rodent) %>% 
           filter(Rodent == 1)
known_rodents <- select(rodents, New.Code) %>% 
                 filter(New.Code != 'DX', New.Code != 'NX',
                        New.Code != 'OX', New.Code != 'PX',
                        New.Code != 'RX', New.Code != 'SX', New.Code != 'UR')
colnames(known_rodents)[1] <- "species"

surveys <- surveys %>% 
           filter(period > 0 & period <= 436)

portal <- semi_join(surveys, known_rodents, by = "species")

# 