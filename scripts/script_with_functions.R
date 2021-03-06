### Cleaned data
# EKB
# March 7, 2016

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

find_reg_persist <- function(data){
  # create column for relative regional persistence (years overall)
  dat <- select(data, species, year) %>% 
         group_by(species, year) %>% 
         summarise(year_count = n())
  dat1 <- select(dat, species) %>% 
          group_by(species) %>% 
          summarise(years = n())
  total_years <- length(unique(data$year))
  dat1 <- mutate(dat1, rel_reg_persist = years/total_years)
  return(dat1)
}

find_reg_occup <- function(data){
  # create column for relative regional occupancy (sites overall)
  dat <- select(data, species, siteID) %>% 
         group_by(species, siteID) %>% 
         summarise(site_count = n())
  dat1 <- select(dat, species) %>% 
          group_by(species) %>% 
          summarise(sites = n())
  total_sites <- length(unique(data$siteID))
  dat1 <- mutate(dat1, rel_reg_occup = sites/total_sites)
  return(dat1)
}

join_regionals <- function(data){
  # combine the relative regional persistance and regional occupancy into one data frame
  dat <- left_join(x = find_reg_persist(data), y = find_reg_occup(data), by = "species")
  return(dat)
}

rel_local_occup <- function(data){
  # create a column for average relative local occupancy (site by year)
  dat <- select(data, year, species, siteID) %>% 
         group_by(year, species) %>% 
         summarise(site_by_year = n_distinct(siteID))
  total_sites <- length(unique(data$siteID))
  dat <- mutate(dat, rel_sites_by_year = site_by_year/total_sites)
  dat1 <- dat %>% group_by(species) %>% 
          summarise(mean_local_occup = mean(rel_sites_by_year))
}


rel_local_persist <- function(data){
  # create a column for average relative local persistance (year by site)
  dat <- select(data, year, species, siteID) %>% 
         group_by(siteID, species) %>% 
         summarise(years_by_site = n_distinct(year))
  total_years <- length(unique(data$year))
  dat <- mutate(dat, rel_years_by_site = years_by_site/total_years)
  dat1 <- dat %>% group_by(species) %>% 
          summarize(mean_local_persist = mean(rel_years_by_site))
  return(dat1)
}

join_locals <- function(data){
  # combine average relative local persistance and occupancy into one data frame
  dat <- left_join(x = rel_local_persist(data), y = rel_local_occup(data), by = "species")
  return(dat)
}

combine_all <- function(data){
  # combine regional and local persistance and occupancy into one data frame
  dat <- left_join(x = join_regionals(data), y = join_locals(data), by = "species")
  return(dat)
}

all_together <- function(data){
  # run all functions together to get the full output
  dat <- select_data(data) # can use pipes here
  dat <- group_data(dat)
  dat <- add_siteID(dat)
  dat <- combine_all(dat)
  return(dat)
}

# GRAPHING functions

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")

graph_local <- function(data){
  # graph of mean local persistance as a function of mean local occupancy
  plot <- ggplot(data, aes(x = mean_local_occup, y = mean_local_persist)) +
    geom_point(aes(color = LTER), size = 6) +
    xlab("Mean Local Occupancy (Spatial)") + 
    ylab("Mean Local Persistence (Temporal)") +
    scale_color_manual(values = cbbPalette) +
    theme(panel.background = element_blank(), 
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"), 
          panel.grid.major = element_line(colour = "light gray"),
          axis.title = element_text(size = 28, face = "bold"),
          axis.text = element_text(size = 20),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 20))
  return(plot)
}

graph_regional <- function(data){
  # graph of mean local persistance as a function of mean local occupancy
  plot <- ggplot(data, aes(x = rel_reg_occup, y = rel_reg_persist)) +
    geom_point(aes(color = LTER), size = 6) +
    xlab("Regional Occupancy (Spatial)") + 
    ylab("Regional Persistence (Temporal)") +
    scale_color_manual(values = cbbPalette) +
    theme(panel.background = element_blank(), 
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"), 
          panel.grid.major = element_line(colour = "light gray"),
          axis.title = element_text(size = 28, face = "bold"),
          axis.text = element_text(size = 20),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 20)) 
  return(plot)
}

#####################
# APPLY FUNCTIONS TO DATASETS

# all functions for all datasets

jor_data <- all_together(jor_data) # I can't figure out how to loop through these
sev_data <- all_together(sev_data)
hja_data <- all_together(hja_data)
sgs_data <- all_together(sgs_data)

# add a dataset column
jor_data$LTER <- "Jornada Basin"    # is there a way to do this with a loop?
hja_data$LTER <- "H.J. Andrews"
sev_data$LTER <- "Sevilleta"
sgs_data$LTER <- "Shortgrass Steppe"

# combine all datasets into one
all_data <- bind_rows(jor_data, sev_data, hja_data, sgs_data)

# graph all data
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")

graph_regional(all_data)
ggsave(file = "all_regional.png", width = 13, height = 10.5)
graph_local(all_data)
ggsave(file = "all_local.png", width = 13, height = 10.5)

###################################################################
# WORKING AREA


#all_output <- all_data %>%
#  group_by(dataset) %>%
#  do(all_together)

#datasets = list(jor_data, sev_data, hja_data, sgs_data)
#output = lapply(datasets, all_together)

#output = data.frame()
#for (dataset in datasets){
  #dataset_output <- all_together(dataset)
  #output = rbind(output, dataset_output)
#}

#hja <- select_data(hja_data) %>% add_siteID()
ggthemes_data$colorblind  <- ggthemes_data$colorblind[-1]
assignInNamespace("ggthemes_data", ggthemes_data, ns="ggthemes")

all_data <- all_data %>% mutate(ymin = mean_local_persist - sd_local_persist, 
                                ymax = mean_local_persist + sd_local_persist,
                                xmin = mean_local_occup - sd_local_occup,
                                xmax = mean_local_occup + sd_local_occup)

ggplot(all_data, aes(x = mean_local_occup, y = mean_local_persist)) +
  geom_point(aes(color = LTER), size = 6) +
  geom_errorbar(aes(ymax = ymax, ymin = ymin))+
  geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
  xlab("Local Occupancy") + 
  ylab("Local Persistence") +
  scale_color_manual(values = cbbPalette) +
  theme(panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"),
        axis.title = element_text(size = 28, face = "bold"),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24)) 
dev.copy(png, "plots/local_sd_error.png", width = 800, height = 550)
dev.off()
