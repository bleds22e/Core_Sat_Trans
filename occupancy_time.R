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

get_occupancy <- function(data){
  # combine all relevant functions to create an occupancy dataframe
  dat <- select_data(data) %>% 
    group_data() %>% 
    add_siteID() %>% 
    occupancy()
  return(dat)
}

plot_occupancy <- function(data){
  # plotting occupancy against time with individual lines for species
  plot <- ggplot(data, aes(x = year, y = occupancy, color = species)) +
          geom_point(size = 2) +
          geom_line(size = 1) +
          theme(panel.background = element_blank(), 
                axis.line = element_line(colour = "black"), 
                panel.grid.major = element_line(colour = "light gray"))
  return(plot)
}

#####################
# OCCUPANCIES

hja_occupancy <- get_occupancy(hja_data)
jor_occupancy <- get_occupancy(jor_data)
sev_occupancy <- get_occupancy(sev_data)
sgs_occupancy <- get_occupancy(sgs_data)

#####################
# OCCUPANCY PLOTS

hja_occ_plot <- plot_occupancy(hja_occupancy)
ggsave(file = "occupancy_hja.png", width = 7, height = 6)

jor_occ_plot <- plot_occupancy(jor_occupancy)
ggsave(file = "occupancy_jor.png", width = 7, height = 6)

sev_occ_plot <- plot_occupancy(sev_occupancy)
ggsave(file = "occupancy_sev.png", width = 7, height = 6)

sgs_occ_plot <- plot_occupancy(sgs_occupancy)
ggsave(file = "occupancy_sgs.png", width = 7, height = 6)