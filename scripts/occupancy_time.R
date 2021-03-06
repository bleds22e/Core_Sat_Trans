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

####################
# FUNCTIONS

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

calculate_cv <- function(occupancy_data){
  # calculates the mean and standard deviation for occupancy over time
  # and adds a new column of the coefficients of variation
  cv <- occupancy_data %>% 
        group_by(species) %>% 
        summarise(mean = mean(occupancy), sd = sd(occupancy)) %>% 
        mutate(cv = sd/mean)
  return(cv)
}

plot_cv_occupancy <- function(cv_data){
  # plot of the coefficients of variation against mean occupancy
  cv_plot <- ggplot(cv_data, aes(x = mean, y = cv)) +
                geom_point(size = 2) +
                xlab("Mean Occupancy") +
                ylab("Coefficients of Variation") + 
                theme(panel.background = element_blank(), 
                    axis.line = element_line(colour = "black"), 
                    panel.grid.major = element_blank())
  return(cv_plot)
}

plot_sd_occupancy <- function(cv_data){
  # plot of the coefficients of variation against mean occupancy
  sd_plot <- ggplot(cv_data, aes(x = mean, y = sd)) +
    geom_point(size = 2) +
    xlab("Mean Occupancy") +
    ylab("Coefficients of Variation") + 
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank())
  return(sd_plot)
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

#####################
# COEFFICIENTS OF VARIATION

hja_cv <- calculate_cv(hja_occupancy)
jor_cv <- calculate_cv(jor_occupancy)
sev_cv <- calculate_cv(sev_occupancy)
sgs_cv <- calculate_cv(sgs_occupancy)

######################
# CV PLOTS

hja_cv_plot <- plot_cv_occupancy(hja_cv)
jor_cv_plot <- plot_cv_occupancy(jor_cv)
sev_cv_plot <- plot_cv_occupancy(sev_cv)
sgs_cv_plot <- plot_cv_occupancy(sgs_cv)

hja_cv_plot
jor_cv_plot
sev_cv_plot
sgs_cv_plot

hja_sd_plot <- plot_sd_occupancy(hja_cv)
jor_sd_plot <- plot_sd_occupancy(jor_cv)
sev_sd_plot <- plot_sd_occupancy(sev_cv)
sgs_sd_plot <- plot_sd_occupancy(sgs_cv)

hja_sd_plot
jor_sd_plot
sev_sd_plot
sgs_sd_plot

#########################
# WORKING

# add a dataset column
jor_cv$LTER <- "Jornada Basin"  
hja_cv$LTER <- "H.J. Andrews"
sev_cv$LTER <- "Sevilleta"
sgs_cv$LTER <- "Shortgrass Steppe"

all_cv_data <- rbind(hja_cv, jor_cv, sev_cv, sgs_cv)

ggplot(all_cv_data, aes(x = mean, y = cv, color = LTER)) +
    geom_point(size = 4) +
    geom_smooth(se = FALSE) +
    xlab("Mean Occupancy") +
    ylab("Coefficients of Variation") + 
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          panel.grid.major = element_line(colour = "light gray"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12)) 

cv_plot <- ggplot(all_cv_data, aes(x = mean, y = cv, color = LTER)) +
  geom_point(size = 4) +
  stat_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
  xlab("Mean Occupancy") +
  ylab("Coefficients of Variation") + 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12)) 
ggsave(cv_plot, file = "all_cv_plot.png", width = 7, height = 6)
