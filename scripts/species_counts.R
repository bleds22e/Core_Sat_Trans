# Ellen Bledsoe
# starting Jan 19, 2016

# this code is to determine the number of species that fall into four different categories:
  # core(spatial) - core(temporal)
  # core(spatial) - transient
  # core(temporal) - satellite
  # transient - satellite

# libraries

library(dplyr)
library(tidyr)
library(ggplot2)

# helpful functions

core_or_sat <- function(rel_sites){
  # function determining whether species is core or satellite (spatial)
  if (rel_sites > 0.5){
    core_sat = 'core'
  } else {
    core_sat = 'sat'
  }
  return(core_sat)
}

core_or_trans <- function(rel_years){
  # function determining whether species is core or transient (temporal)
  if (rel_years > 0.5){
    core_trans = 'core'
  } else {
    core_trans = 'trans'
  }
  return(core_trans)
}

category_3x3 <- function(rel_sites, rel_years){
  if (rel_sites >= (2/3) & rel_years >= (2/3)){
    category = '1'
  } else if (rel_sites >= (2/3) & rel_years < (2/3) & rel_years > (1/3)){
    category = '2'
  } else if (rel_sites >= (2/3) & rel_years <= (1/3)){
    category = '3'
  } else if (rel_sites < (2/3) & rel_sites > (1/3) & rel_years >= (2/3)){
    category = '4'
  } else if (rel_sites < (2/3) & rel_sites > (1/3) & rel_years < (2/3) & rel_years > (1/3)){
    category = '5'
  } else if (rel_sites < (2/3) & rel_sites > (1/3) & rel_years <= (1/3)){
    category = '6'
  } else if (rel_sites <= (1/3) & rel_years >= (2/3)){
    category = '7'
  } else if (rel_sites <= (1/3) & rel_years < (2/3) & rel_years > (1/3)){
    category = '8'
  } else {
    category = '9'
  }
  return(category)
}

# read in files
jornada <- read.csv("data/Rodents/jornada_rodents.csv", header = TRUE, na.strings = ".")
head(jornada)

sevilleta <- read.csv("data/Rodents/sevilleta_sm_mrc.txt", header = TRUE)
head(sevilleta)

hj_andrews <- read.csv("data/Rodents/hj_andrews.csv", header = TRUE)
head(hj_andrews)

shortgrass_steppe <- read.table("data/Rodents/SGS_LTER_smammals.txt", header = TRUE, sep = "\t", na.strings = ".")
head(shortgrass_steppe)

# unique trap sessions
jor <- select(jornada, year, season, habitat) %>% 
       group_by(year, season, habitat) %>% 
       arrange(year, season, habitat)
jor_trap_sessions <- unique(jor)
jor_trap_sessions

sev <- select(sevilleta, year, season, location) %>% 
       group_by(year, season, location) %>% 
       arrange(year, season, location)
sev_trap_sessions <- unique(sev)
sev_trap_sessions

hja <- select(hj_andrews, YEAR, REPBLK, SUBPLOT) %>% 
       group_by(YEAR, REPBLK, SUBPLOT) %>% 
       arrange(YEAR, REPBLK, SUBPLOT) 
hja_trap_sessions <- unique(hja)
hja_trap_sessions

sgs <- select(shortgrass_steppe, YEAR, VEG, WEB) %>% 
  group_by(YEAR, VEG, WEB) %>% 
  arrange(YEAR, VEG, WEB) 
sgs_trap_sessions <- unique(sgs)
sgs_trap_sessions

### extracting the requisite data

# jornada
jor_data <- select(jornada, year, season, habitat, web, spp, recap) %>% 
            filter(recap != "Y", spp != "DIPO1", spp != "PERO1", spp != "NA") %>% 
            group_by(spp, year, season, habitat, web) %>% 
            arrange(spp, year, season, habitat, web) %>% 
            summarise(count = n())
# give each site a unique id
jor_data$site_ID <- paste(jor_data$habitat, jor_data$web, sep = "_")

### number of years in which species occur

jor_years <- select(jor_data, spp, year) %>% 
             group_by(spp, year) %>% 
             summarise(year_count = n())
jor_years_count <- select(jor_years, spp) %>% 
                   group_by(spp) %>% 
                   summarise(years = n())
# total years
total_years <- length(unique(jor_data$year))
# histogram of years
hist(jor_years_count$years)

# add column for relative number of years present
jor_years_count <- mutate(jor_years_count, rel_years = years/total_years)
jor_years_count

### number of locations in which species occur

jor_locs <- select(jor_data, spp, site_ID) %>% 
            group_by(spp, site_ID) %>% 
            summarise(site_count = n())
jor_locs_count <- select(jor_locs, spp) %>% 
                  group_by(spp) %>% 
                  summarise(sites = n())
# total locations
total_sites <- length(unique(jor_data$site_ID))

# add column for relative number of sites present
jor_locs_count <- mutate(jor_locs_count, rel_sites = sites/total_sites)

### combine all relevant data into one dataframe

jor_rel_data <- left_join(jor_years_count, jor_locs_count, by = "spp")

jor_rel_years_list <- as.list(jor_rel_data$rel_years)
jor_rel_sites_list <- as.list(jor_rel_data$rel_sites)
jor_rel_data$core_sat <- rapply(jor_rel_sites_list, core_or_sat)
jor_rel_data$core_trans <- rapply(jor_rel_years_list, core_or_trans)
jor_rel_data$category <- paste(jor_rel_data$core_sat, jor_rel_data$core_trans, sep = "_")
jor_rel_data$cat_3x3 <- mapply(category_3x3, jor_rel_data$rel_years, jor_rel_data$rel_sites)

jor_rel_data

### SEVILLETA ###

# sevilleta
sev_data <- select(sevilleta, year, season, location, web, species, recap) %>% 
            filter(recap != "y", species != "pgsp", species != "dipo", 
               species != "nesp", species != "onsp", species != "pesp", 
               species != "resp", species != "na", species != "pmsp", species != "spsp",
               year > 1991 & year < 2007, 
               location != "two22", location != "blugrama", location != "savanna") %>% 
            group_by(species, year, season, location, web) %>% 
            arrange(species, year, season, location, web) %>% 
            summarise(count = n())
# unique location ids
sev_data$site_ID <- paste(sev_data$location, sev_data$web, sep = "_")

sev_years <- select(sev_data, species, year) %>% 
  group_by(species, year) %>% 
  summarise(year_count = n())
sev_years_count <- select(sev_years, species) %>% 
  group_by(species) %>% 
  summarise(years = n())
# total years
total_years_sev <- length(unique(sev_data$year))
# histogram of years
hist(sev_years_count$years)

# add column for relative number of years present
sev_years_count <- mutate(sev_years_count, rel_years = years/total_years_sev)
sev_years_count

### number of locations in which species occur

sev_locs <- select(sev_data, species, site_ID) %>% 
  group_by(species, site_ID) %>%
  summarise(site_count = n())
sev_locs_count <- select(sev_locs, species) %>% 
  group_by(species) %>% 
  summarise(sites = n())
# total locations
total_sites_sev <- length(unique(sev_data$site_ID))
total_sites_sev
# histogram of site frequency
hist(sev_locs_count$sites)

# add column for relative number of sites present
sev_locs_count <- mutate(sev_locs_count, rel_sites = sites/total_sites_sev)

# combine all relevant data into one dataframe

sev_rel_data <- left_join(sev_years_count, sev_locs_count, by = "species")

rel_years_list <- as.list(sev_rel_data$rel_years)
rel_sites_list <- as.list(sev_rel_data$rel_sites)
sev_rel_data$core_sat <- rapply(rel_sites_list, core_or_sat)
sev_rel_data$core_trans <- rapply(rel_years_list, core_or_trans)
sev_rel_data$category <- paste(sev_rel_data$core_sat, sev_rel_data$core_trans, sep = "_")
sev_rel_data$cat_3x3 <- mapply(category_3x3, sev_rel_data$rel_years, sev_rel_data$rel_sites)

sev_rel_data

### different time frame and sites
sev_data1 <- select(sevilleta, year, season, location, web, species, recap) %>% 
             filter(recap != "y", species != "pgsp", species != "dipo", 
                    species != "nesp", species != "onsp", species != "pesp", 
                    species != "resp", species != "na", species != "pmsp", species != "spsp", 
                    year < 2009, location != "two22", location != "blugrama", 
                    location != "savanna", location != "goatdraw", location != "rsgrass") %>% 
             group_by(species, year, season, location, web) %>% 
             arrange(species, year, season, location, web) %>% 
             summarise(count = n())
# unique location ids
sev_data1$site_ID <- paste(sev_data1$location, sev_data1$web, sep = "_")
sev_data1

sev_years1 <- select(sev_data1, species, year) %>% 
  group_by(species, year) %>% 
  summarise(year_count = n())
sev_years_count1 <- select(sev_years1, species) %>% 
  group_by(species) %>% 
  summarise(years = n())
# total years
total_years_sev1 <- length(unique(sev_data1$year))
# histogram of years
hist(sev_years_count1$years)

# add column for relative number of years present
sev_years_count1 <- mutate(sev_years_count1, rel_years = years/total_years_sev1)
sev_years_count1

# number of locations in which species occur

sev_locs1 <- select(sev_data1, species, site_ID) %>% 
  group_by(species, site_ID) %>%
  summarise(site_count = n())
sev_locs_count1 <- select(sev_locs1, species) %>% 
  group_by(species) %>% 
  summarise(sites = n())
# total locations
total_sites_sev1 <- length(unique(sev_data1$site_ID))
total_sites_sev1
# histogram of site frequency
hist(sev_locs_count1$sites)

# add column for relative number of sites present
sev_locs_count1 <- mutate(sev_locs_count1, rel_sites = sites/total_sites_sev1)

# combine all relevant data into one dataframe

sev_rel_data1 <- left_join(sev_years_count1, sev_locs_count1, by = "species")

rel_years_list <- as.list(sev_rel_data1$rel_years)
rel_sites_list <- as.list(sev_rel_data1$rel_sites)
sev_rel_data1$core_sat <- rapply(rel_sites_list, core_or_sat)
sev_rel_data1$core_trans <- rapply(rel_years_list, core_or_trans)
sev_rel_data1$category <- paste(sev_rel_data1$core_sat, sev_rel_data1$core_trans, sep = "_")
sev_rel_data1$cat_3x3 <- mapply(category_3x3, sev_rel_data1$rel_years, sev_rel_data1$rel_sites)

sev_rel_data1

### HJ Andrews ###

# hj andrews
hja_data <- select(hj_andrews, YEAR, REPBLK, SUBPLOT, SPECIES) %>% 
  filter(SPECIES != 'UNKN') %>% 
  group_by(SPECIES, YEAR, REPBLK, SUBPLOT) %>% 
  arrange(SPECIES, YEAR, REPBLK, SUBPLOT) %>% 
  summarise(count = n())
# unique location ids
hja_data$site_ID <- paste(hja_data$REPBLK, hja_data$SUBPLOT, sep = "_")

hja_years <- select(hja_data, SPECIES, YEAR) %>% 
  group_by(SPECIES, YEAR) %>% 
  summarise(year_count = n())
hja_years_count <- select(hja_years, SPECIES) %>% 
  group_by(SPECIES) %>% 
  summarise(years = n())
# total years
total_years_hja <- length(unique(hja_data$YEAR))
# histogram of years
hist(hja_years_count$years)

# add column for relative number of years present
hja_years_count <- mutate(hja_years_count, rel_years = years/total_years_hja)
hja_years_count

### number of locations in which species occur

hja_locs <- select(hja_data, SPECIES, site_ID) %>% 
  group_by(SPECIES, site_ID) %>%
  summarise(site_count = n())
hja_locs_count <- select(hja_locs, SPECIES) %>% 
  group_by(SPECIES) %>% 
  summarise(sites = n())
# total locations
total_sites_hja <- length(unique(hja_data$site_ID))
total_sites_hja
# histogram of site frequency
hist(hja_locs_count$sites)

# add column for relative number of sites present
hja_locs_count <- mutate(hja_locs_count, rel_sites = sites/total_sites_hja)

# combine all relevant data into one dataframe

hja_rel_data <- left_join(hja_years_count, hja_locs_count, by = "SPECIES")

hja_rel_years_list <- as.list(hja_rel_data$rel_years)
hja_rel_sites_list <- as.list(hja_rel_data$rel_sites)
hja_rel_data$core_sat <- rapply(hja_rel_sites_list, core_or_sat)
hja_rel_data$core_trans <- rapply(hja_rel_years_list, core_or_trans)
hja_rel_data$category_2x2 <- paste(hja_rel_data$core_sat, hja_rel_data$core_trans, sep = "_")
hja_rel_data$cat_3x3 <- mapply(category_3x3, hja_rel_data$rel_years, hja_rel_data$rel_sites)
hja_rel_data

### Shortgrass Steppe ###

sgs_species <- select(shortgrass_steppe, SPP)
sgs_species <- unique(sgs_species)

# sgs
sgs_data <- select(shortgrass_steppe, YEAR, VEG, WEB, SPP) %>% 
  filter(SPP != 'NA') %>% 
  group_by(SPP, YEAR, VEG, WEB) %>% 
  arrange(SPP, YEAR, VEG, WEB) %>% 
  summarise(count = n())
# unique location ids
sgs_data$site_ID <- paste(sgs_data$VEG, sgs_data$WEB, sep = "_")

sgs_years <- select(sgs_data, SPP, YEAR) %>% 
  group_by(SPP, YEAR) %>% 
  summarise(year_count = n())
sgs_years_count <- select(sgs_years, SPP) %>% 
  group_by(SPP) %>% 
  summarise(years = n())
# total years
total_years_sgs <- length(unique(sgs_data$YEAR))
# histogram of years
hist(sgs_years_count$years)

# add column for relative number of years present
sgs_years_count <- mutate(sgs_years_count, rel_years = years/total_years_sgs)
sgs_years_count

### number of locations in which species occur

sgs_locs <- select(sgs_data, SPP, site_ID) %>% 
  group_by(SPP, site_ID) %>%
  summarise(site_count = n())
sgs_locs_count <- select(sgs_locs, SPP) %>% 
  group_by(SPP) %>% 
  summarise(sites = n())
# total locations
total_sites_sgs <- length(unique(sgs_data$site_ID))
total_sites_sgs
# histogram of site frequency
hist(sgs_locs_count$sites)

# add column for relative number of sites present
sgs_locs_count <- mutate(sgs_locs_count, rel_sites = sites/total_sites_sgs)

# combine all relevant data into one dataframe

sgs_rel_data <- left_join(sgs_years_count, sgs_locs_count, by = "SPP")

sgs_rel_years_list <- as.list(sgs_rel_data$rel_years)
sgs_rel_sites_list <- as.list(sgs_rel_data$rel_sites)
sgs_rel_data$core_sat <- rapply(sgs_rel_sites_list, core_or_sat)
sgs_rel_data$core_trans <- rapply(sgs_rel_years_list, core_or_trans)
sgs_rel_data$category_2x2 <- paste(sgs_rel_data$core_sat, sgs_rel_data$core_trans, sep = "_")
sgs_rel_data$cat_3x3 <- mapply(category_3x3, sgs_rel_data$rel_years, sgs_rel_data$rel_sites)
sgs_rel_data

###

hja_rel_data
sev_rel_data1
sev_rel_data
sgs_rel_data
jor_rel_data

### combine and graph
colnames(jor_rel_data)[1] <- "species"
colnames(hja_rel_data)[1] <- "species"
colnames(sev_rel_data)[1] <- "species"
colnames(sev_rel_data1)[1] <- "species"
colnames(sgs_rel_data)[1] <- "species"

colnames(jor_rel_data)[8] <- "cat_2x2"
colnames(hja_rel_data)[8] <- "cat_2x2"
colnames(sev_rel_data)[8] <- "cat_2x2"
colnames(sev_rel_data1)[8] <- "cat_2x2"
colnames(sgs_rel_data)[8] <- "cat_2x2"

jor_rel_data$dataset <- "jor"
hja_rel_data$dataset <- "hja"
sev_rel_data1$dataset <- "sev1"
sev_rel_data$dataset <- "sev"
sgs_rel_data$dataset <- "sgs"

all_data <- bind_rows(hja_rel_data, sev_rel_data1, sgs_rel_data, jor_rel_data)
head(all_data)

plot(all_data$rel_years ~ all_data$rel_sites)

ggplot(all_data, aes(x = rel_sites, y = rel_years)) +
  geom_point(aes(color = dataset), size = 3) +
  xlab("Relative Number of Sites Occupied") + 
  ylab("Relative Number of Years Present") +
  ggtitle("Rel. Years (overall) ~ Rel. Sites (overall)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "all_broad_yr_site.png", width = 6, height = 5)

### Try relative number of sites per year---keep years intact 
  # overall years, parsing space by time

## jornada

jor_data1 <- select(jornada, year, season, habitat, web, spp, recap) %>% 
             filter(recap != "Y", spp != "DIPO1", spp != "PERO1", spp != "NA") %>% 
             group_by(spp, year, season, habitat, web) %>% 
             arrange(spp, year, season, habitat, web)

# give each site a unique id
jor_data1$site_ID <- paste(jor_data1$habitat, jor_data1$web, sep = "_")

jor_sites_by_year <- select(jor_data1, year, spp, site_ID) %>% 
                 group_by(year, spp) %>% 
                 summarise(site_by_year = n_distinct(site_ID))

# total years
jor_total_sites <- length(unique(jor_data1$site_ID))

# add column for relative number of sites per year
jor_rel_sites_by_year <- mutate(jor_sites_by_year, rel_sites_by_year = site_by_year/jor_total_sites)
jor_rel_sites_by_year

# average relative sites per year by species
jor_sites_year_species <- jor_rel_sites_by_year %>% 
                          group_by(spp) %>% 
                          summarize(mean_rel_sites = mean(rel_sites_by_year))
jor_sites_year_species
colnames(jor_sites_year_species)[1] <- "species"
colnames(jor_rel_data)[1] <- "species"
jor_rel_data_plus <- inner_join(jor_rel_data, jor_sites_year_species, by = "species")
jor_rel_data_plus
hist(jor_rel_data_plus$rel_sites)
hist(jor_rel_data_plus$mean_rel_sites)
plot(jor_rel_data_plus$rel_years ~ jor_rel_data_plus$mean_rel_sites)

## sevilleta

sev_data1_1 <- select(sevilleta, year, season, location, web, species, recap) %>% 
  filter(recap != "y", species != "pgsp", species != "dipo", 
         species != "nesp", species != "onsp", species != "pesp", 
         species != "resp", species != "na", species != "pmsp", species != "spsp", 
         year < 2009, location != "two22", location != "blugrama", 
         location != "savanna", location != "goatdraw") %>% 
  group_by(species, year, season, location, web) %>% 
  arrange(species, year, season, location, web) 
  
# unique location ids
sev_data1_1$site_ID <- paste(sev_data1_1$location, sev_data1_1$web, sep = "_")
sev_data1_1

sev_sites_by_year <- select(sev_data1_1, year, species, site_ID) %>% 
  group_by(year, species) %>% 
  summarise(site_by_year = n_distinct(site_ID))

# total years
sev_total_sites <- length(unique(sev_data1_1$site_ID))

# add column for relative number of sites per year
sev_rel_sites_by_year <- mutate(sev_sites_by_year, rel_sites_by_year = site_by_year/sev_total_sites)
sev_rel_sites_by_year

# average relative sites per year by species
sev_sites_year_species <- sev_rel_sites_by_year %>% 
  group_by(species) %>% 
  summarize(mean_rel_sites = mean(rel_sites_by_year))
sev_sites_year_species
colnames(sev_sites_year_species)[1] <- "species"
sev_rel_data_plus <- inner_join(sev_rel_data, sev_sites_year_species, by = "species")
sev_rel_data_plus
plot(sev_rel_data_plus$rel_years ~ sev_rel_data_plus$mean_rel_sites)

## hj andrews

hja_data1 <- select(hj_andrews, YEAR, REPBLK, SUBPLOT, SPECIES) %>% 
  filter(SPECIES != 'UNKN') %>% 
  group_by(SPECIES, YEAR, REPBLK, SUBPLOT) %>% 
  arrange(SPECIES, YEAR, REPBLK, SUBPLOT) 
 
# unique location ids
hja_data1$site_ID <- paste(hja_data1$REPBLK, hja_data1$SUBPLOT, sep = "_")

hja_sites_by_year <- select(hja_data1, YEAR, SPECIES, site_ID) %>% 
  group_by(YEAR, SPECIES) %>% 
  summarise(site_by_year = n_distinct(site_ID))

# total years
hja_total_sites <- length(unique(hja_data1$site_ID))

# add column for relative number of sites per year
hja_rel_sites_by_year <- mutate(hja_sites_by_year, rel_sites_by_year = site_by_year/hja_total_sites)
hja_rel_sites_by_year

# average relative sites per year by species
hja_sites_year_species <- hja_rel_sites_by_year %>% 
                          group_by(SPECIES) %>% 
                          summarize(mean_rel_sites = mean(rel_sites_by_year))
hja_sites_year_species
colnames(hja_sites_year_species)[1] <- "species"
hja_rel_data_plus <- inner_join(hja_rel_data, hja_sites_year_species, by = "species")
hja_rel_data_plus
plot(hja_rel_data_plus$rel_years ~ hja_rel_data_plus$mean_rel_sites)

## shortgrass steppe

sgs_data1 <- select(shortgrass_steppe, YEAR, VEG, WEB, SPP) %>% 
            filter(SPP != 'NA') %>% 
            group_by(SPP, YEAR, VEG, WEB) %>% 
            arrange(SPP, YEAR, VEG, WEB)

# unique location ids
sgs_data1$site_ID <- paste(sgs_data1$VEG, sgs_data1$WEB, sep = "_")

sgs_sites_by_year <- select(sgs_data1, YEAR, SPP, site_ID) %>% 
  group_by(YEAR, SPP) %>% 
  summarise(site_by_year = n_distinct(site_ID))

# total years
sgs_total_sites <- length(unique(sgs_data1$site_ID))

# add column for relative number of sites per year
sgs_rel_sites_by_year <- mutate(sgs_sites_by_year, rel_sites_by_year = site_by_year/sgs_total_sites)
sgs_rel_sites_by_year

# average relative sites per year by species
sgs_sites_year_species <- sgs_rel_sites_by_year %>% 
  group_by(SPP) %>% 
  summarize(mean_rel_sites = mean(rel_sites_by_year))
sgs_sites_year_species
colnames(sgs_sites_year_species)[1] <- "species"
sgs_rel_data_plus <- inner_join(sgs_rel_data, sgs_sites_year_species, by = "species")
sgs_rel_data_plus
plot(sgs_rel_data_plus$rel_years ~ sgs_rel_data_plus$mean_rel_sites)

## combine everything

all_data_plus <- bind_rows(hja_rel_data_plus, sev_rel_data_plus, sgs_rel_data_plus, jor_rel_data_plus)
head(all_data_plus)

plot(all_data_plus$rel_years ~ all_data_plus$mean_rel_sites)

ggplot(all_data_plus, aes(x = mean_rel_sites, y = rel_years)) +
  geom_point(aes(color = dataset), size = 3) +
  xlab("Relative Number of Sites Occupied Per Year") + 
  ylab("Relative Number of Years Present") +
  ggtitle("Rel. Years (overall) ~ Rel. Sites (per year)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "all_mean_rel_site.png", width = 6, height = 5)

### Try relative number of years per site---keep sites intact 
  # overall sites, parsing space by time

## jornada

jor_data1 <- select(jornada, year, season, habitat, web, spp, recap) %>% 
             filter(recap != "Y", spp != "DIPO1", spp != "PERO1", spp != "NA") %>% 
             group_by(spp, year, season, habitat, web) %>% 
             arrange(spp, year, season, habitat, web)

# give each site a unique id
jor_data1$site_ID <- paste(jor_data1$habitat, jor_data1$web, sep = "_")

jor_years_by_site <- select(jor_data1, year, spp, site_ID) %>% 
                     group_by(site_ID, spp) %>% 
                     summarise(years_by_site = n_distinct(year))

# total years
jor_total_years <- length(unique(jor_data1$year))
jor_total_years

# add column for relative number of sites per year
jor_rel_years_by_site <- mutate(jor_years_by_site, rel_years_by_site = years_by_site/jor_total_years)
jor_rel_years_by_site

# average relative sites per year by species
jor_years_site_species <- jor_rel_years_by_site %>% 
                          group_by(spp) %>% 
                          summarize(mean_rel_years = mean(rel_years_by_site))
jor_years_site_species
colnames(jor_years_site_species)[1] <- "species"
jor_rel_data_3 <- inner_join(jor_rel_data_plus, jor_years_site_species, by = "species")
jor_rel_data_3

plot(jor_rel_data_3$mean_rel_years ~ jor_rel_data_3$rel_sites)
plot(jor_rel_data_3$mean_rel_years ~ jor_rel_data_3$mean_rel_sites)

## sevilleta

sev_data1_1 <- select(sevilleta, year, season, location, web, species, recap) %>% 
  filter(recap != "y", species != "pgsp", species != "dipo", 
         species != "nesp", species != "onsp", species != "pesp", 
         species != "resp", species != "na", species != "pmsp", species != "spsp", 
         year < 2009, location != "two22", location != "blugrama", 
         location != "savanna", location != "goatdraw") %>% 
  group_by(species, year, season, location, web) %>% 
  arrange(species, year, season, location, web) 

# unique location ids
sev_data1_1$site_ID <- paste(sev_data1_1$location, sev_data1_1$web, sep = "_")
sev_data1_1

sev_years_by_site <- select(sev_data1_1, year, species, site_ID) %>% 
  group_by(site_ID, species) %>% 
  summarise(years_by_site = n_distinct(year))

# total years
sev_total_years <- length(unique(sev_data1_1$year))
sev_total_years

# add column for relative number of sites per year
sev_rel_years_by_site <- mutate(sev_years_by_site, rel_years_by_site = years_by_site/sev_total_years)
sev_rel_years_by_site

# average relative sites per year by species
sev_years_site_species <- sev_rel_years_by_site %>% 
                          group_by(species) %>% 
                          summarize(mean_rel_years = mean(rel_years_by_site))
sev_years_site_species
colnames(sev_years_site_species)[1] <- "species"
sev_rel_data_3 <- inner_join(sev_rel_data_plus, sev_years_site_species, by = "species")
sev_rel_data_3

plot(sev_rel_data_3$mean_rel_years ~ sev_rel_data_3$rel_sites)
plot(sev_rel_data_3$mean_rel_years ~ sev_rel_data_3$mean_rel_sites)

## hj andrews

hja_data1 <- select(hj_andrews, YEAR, REPBLK, SUBPLOT, SPECIES) %>% 
  filter(SPECIES != 'UNKN') %>% 
  group_by(SPECIES, YEAR, REPBLK, SUBPLOT) %>% 
  arrange(SPECIES, YEAR, REPBLK, SUBPLOT) 

# unique location ids
hja_data1$site_ID <- paste(hja_data1$REPBLK, hja_data1$SUBPLOT, sep = "_")

hja_years_by_site <- select(hja_data1, YEAR, SPECIES, site_ID) %>% 
  group_by(site_ID, SPECIES) %>% 
  summarise(years_by_site = n_distinct(YEAR))

# total years
hja_total_years <- length(unique(hja_data1$YEAR))
hja_total_years

# add column for relative number of sites per year
hja_rel_years_by_site <- mutate(hja_years_by_site, rel_years_by_site = years_by_site/hja_total_years)
hja_rel_years_by_site

# average relative sites per year by species
hja_years_site_species <- hja_rel_years_by_site %>% 
                          group_by(SPECIES) %>% 
                          summarize(mean_rel_years = mean(rel_years_by_site))
hja_years_site_species
colnames(hja_years_site_species)[1] <- "species"
hja_rel_data_3 <- inner_join(hja_rel_data_plus, hja_years_site_species, by = "species")
hja_rel_data_3

plot(hja_rel_data_3$mean_rel_years ~ hja_rel_data_3$rel_sites)
plot(hja_rel_data_3$mean_rel_years ~ hja_rel_data_3$mean_rel_sites)

## shortgrass steppe

sgs_data1 <- select(shortgrass_steppe, YEAR, VEG, WEB, SPP) %>% 
  filter(SPP != 'NA') %>% 
  group_by(SPP, YEAR, VEG, WEB) %>% 
  arrange(SPP, YEAR, VEG, WEB)

# unique location ids
sgs_data1$site_ID <- paste(sgs_data1$VEG, sgs_data1$WEB, sep = "_")

sgs_years_by_site <- select(sgs_data1, YEAR, SPP, site_ID) %>% 
                     group_by(site_ID, SPP) %>% 
                     summarise(years_by_site = n_distinct(YEAR))

# total years
sgs_total_years <- length(unique(sgs_data1$YEAR))
sgs_total_years

# add column for relative number of sites per year
sgs_rel_years_by_site <- mutate(sgs_years_by_site, rel_years_by_site = years_by_site/sgs_total_years)
sgs_rel_years_by_site

# average relative sites per year by species
sgs_years_site_species <- sgs_rel_years_by_site %>% 
                          group_by(SPP) %>% 
                          summarize(mean_rel_years = mean(rel_years_by_site))
sgs_years_site_species
colnames(sgs_years_site_species)[1] <- "species"
sgs_rel_data_3 <- inner_join(sgs_rel_data_plus, sgs_years_site_species, by = "species")
sgs_rel_data_3

plot(sgs_rel_data_3$mean_rel_years ~ sgs_rel_data_3$rel_sites)
plot(sgs_rel_data_3$mean_rel_years ~ sgs_rel_data_3$mean_rel_sites)

## combine everything

all_data_3 <- bind_rows(hja_rel_data_3, sev_rel_data_3, sgs_rel_data_3, jor_rel_data_3)
head(all_data_3)

ggplot(all_data_3, aes(x = rel_sites, y = mean_rel_years)) +
  geom_point(aes(color = dataset), size = 3) +
  xlab("Relative Number of Sites Occupied") + 
  ylab("Relative Number of Years Present Per Site") +
  ggtitle("Rel. Years (per site) ~ Rel. Sites (overall)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "all_mean_rel_year.png", width = 6, height = 5)

ggplot(all_data_3, aes(x = mean_rel_sites, y = mean_rel_years)) +
  geom_point(aes(color = dataset), size = 3) +
  xlab("Relative Number of Sites Occupied Per Year") + 
  ylab("Relative Number of Years Present Per Site") +
  ggtitle("Rel. Years (per site) ~ Rel. Sites (per year)") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "all_mean_rel_year&site.png", width = 6, height = 5)

