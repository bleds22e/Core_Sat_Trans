# Spatiotemporal variation with Sevilleta smammal data
# January 2018
# EKB

library(dplyr)

# read in data
sev_rodents <- read.csv("data/Rodents/sevilleta_sm_mrc.txt", header = TRUE, stringsAsFactors = FALSE)

### Preliminary Data Prep ###

# remove recaps and unknown species
sev_rodents <- sev_rodents %>% filter(recap == "y", species != "pgsp" | species != "dipo" | species != "nesp" | species != "onsp" | species != "pmsp" | species != "resp" | species != "na")

# find number of webs by year and location
unique_webs <- sev_rodents %>% select(year, location, web) %>% unique()

# web counts by year and location
web_counts <- unique_webs %>% group_by(location, web) %>% 
  summarise(count = n(), first_year = min(year), last_year = max(year))

sampling_20plus <- web_counts %>% filter(count >= 20)
sampling_15plus <- web_counts %>% filter(count >= 15)

# select data from locations and webs with 15 or more years of data
data15 <- semi_join(sev_rodents, sampling_15plus, by = c("location", "web"))
data20 <- semi_join(sev_rodents, sampling_20plus, by = c("location", "web"))

### Core-Sporadic Species ###

# by web, by location, by system
# only take same years for all groupings?

data15 <- data15 %>% filter(year >)
data20 <- data20 %>% filter(year < 2010)
