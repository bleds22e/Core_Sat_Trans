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
  filter(SPP != 'NA', SPP != "SOSP")

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
  dat <- select_data(data) %>% group_data() %>% add_siteID() %>% combine_all()
  return(dat)
}

####################
# ABUNDANCES

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
jor_data1$LTER <- "Jornada Basin"  
hja_data1$LTER <- "H.J. Andrews"
sev_data1$LTER <- "Sevilleta"
sgs_data1$LTER <- "Shortgrass Steppe"

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

####################
# GRAPHS

cbbPalette <- c("#009E73", "#e79f00", "#9ad0f3", "#CC79A7")

total_abund <- ggplot(all_abund_data, aes(x = rel_reg_occup, y = rel_reg_persist)) +
  geom_point(aes(color = LTER, size = adj_abund)) +
  ggtitle("Adjusted Total Abundance") +
  xlab("Regional Occupancy") + 
  ylab("Regional Persistence") +
  scale_color_manual(values = cbbPalette) +
  scale_size(range = c(2,12)) +
  theme(panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 28, face = "bold")) 

avg_abund <- ggplot(all_abund_data, aes(x = rel_reg_occup, y = rel_reg_persist)) +
  geom_point(aes(color = LTER, size = ajd_avg_abund)) +
  ggtitle("Adjusted Average Abundance") +
  xlab("Regional Occupancy") + 
  ylab("Regional Persistence") +
  scale_color_manual(values = cbbPalette) +
  scale_size(range = c(2,12)) +
  theme(panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 28, face = "bold"))

require(cowplot)
plot_grid(total_abund, avg_abund, align = 'h')
ggsave(filename = "abundances.png")

#########################
# WORK AREA

# SE error bars

#rel_local_occup <- function(data){
  # create a column for average relative local occupancy (site by year)
  dat <- select(jor_data1, year, species, siteID) %>% 
    group_by(year, species) %>% 
    summarise(site_by_year = n_distinct(siteID))
  total_sites <- n_distinct(data$siteID)
  dat <- mutate(dat, rel_sites_by_year = site_by_year/total_sites)
  dat1 <- dat %>% group_by(species) %>% 
    summarise(mean_local_occup = mean(rel_sites_by_year), 
              SE_local_occup = (sd(rel_sites_by_year)/sqrt(n(rel_sites_by_year))))

#rel_local_persist
  # create a column for average relative local persistance (year by site)
  dat <- select(data, year, species, siteID) %>% 
    group_by(siteID, species) %>% 
    summarise(years_by_site = n_distinct(year))
  total_years <- length(unique(data$year))
  dat <- mutate(dat, rel_years_by_site = years_by_site/total_years)
  dat1 <- dat %>% group_by(species) %>% 
    summarize(mean_local_persist = mean(rel_years_by_site))
  return(dat1)


