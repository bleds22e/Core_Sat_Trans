# Transients by Season
# EKB
# Nov 2018

########################################################################
# This script is different from what's in the `transient_CCA.r` script
# in the following ways:
#       - different NDVI file from PortalData repo
#       - 1988 - 2014
#       - has month and year columns in combined_data
#       - will have seasons added

# LIBRARIES #
library(RCurl)
library(stringr)
library(zoo)
library(lubridate)
library(vegan)
library(raster)
library(tidyverse)

# FUNCTION #

find_reg_persist <- function(data){
  # create column for relative regional persistence (periods overall)
  dat <- dplyr::select(data, species, period) %>% 
    dplyr::group_by(species, period) %>% 
    dplyr::summarise(period_count = n())
  dat1 <- dplyr::select(dat, species) %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(periods = n())
  total_periods <- length(unique(data$period))
  dat1 <- dplyr::mutate(dat1, rel_reg_persist = periods/total_periods)
  return(dat1)
}

# GET DATA #

rdat <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"), na.strings = '')
trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/Landsat_monthly_NDVI.csv"))

#===============================================================

# get rodent abundance datat
abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                            plots = "all", unknowns = F, shape = "flat", 
                            time = "period", min_plots = 10)

# BREAK UP COMMUNITIES BY LDA

# pre-2000, 2000-2009, post-2009 - using only 3 because second isn't very clean and don't have NDVI data before 1992

comm1 <- filter(abund, period > 75 & period < 261) # pre-2000
comm2 <- filter(abund, period >= 261 & period < 381) # 2000-2009
comm3 <- filter(abund, period >= 381) # post-2009

comm1_persist <- filter(comm1, abundance > 0) %>% 
  find_reg_persist()
comm2_persist <- filter(comm2, abundance > 0) %>% 
  find_reg_persist()
comm3_persist <- filter(comm3, abundance > 0) %>% 
  find_reg_persist()

par(mfrow = c(3,1))
hist(comm1_persist$rel_reg_persist)
hist(comm2_persist$rel_reg_persist)
hist(comm3_persist$rel_reg_persist)

#dev.copy(png, "plots/hist_portal_comms.png")
dev.off()

# transients for each community

transients1 <- comm1_persist %>% 
  filter(rel_reg_persist <= 1/3, species != "PB") %>% 
  dplyr::select(species)
transients1 <- unlist(apply(transients1, 1, list), recursive = FALSE)

transients2 <- comm2_persist %>% 
  filter(rel_reg_persist <= 1/3 ) %>% 
  dplyr::select(species)
transients2 <- unlist(apply(transients2, 1, list), recursive = FALSE)

transients3 <- comm3_persist %>% 
  filter(rel_reg_persist <= 1/3) %>% 
  dplyr::select(species)
transients3 <- unlist(apply(transients3, 1, list), recursive = FALSE)


### GET ABUNDANCES ###

# get total abundance by community
total_by_period1 <- comm1 %>% group_by(period) %>% 
  summarise(total = sum(abundance))
total_by_period2 <- comm2 %>% group_by(period) %>% 
  summarise(total = sum(abundance))
total_by_period3 <- comm3 %>% group_by(period) %>% 
  summarise(total = sum(abundance))

# get relative abundances 
rel_abund1 <- filter(comm1, species %in% transients1) %>% 
  group_by(period) %>% 
  summarise(transients = sum(abundance))
rel_abund1 <- inner_join(rel_abund1, total_by_period1)
rel_abund_transients1 <- mutate(rel_abund1, transient_rel_abund = transients/total)

rel_abund2 <- filter(comm2, species %in% transients2) %>% 
  group_by(period) %>% 
  summarise(transients = sum(abundance))
rel_abund2 <- inner_join(rel_abund2, total_by_period2)
rel_abund_transients2 <- mutate(rel_abund2, transient_rel_abund = transients/total)

rel_abund3 <- filter(comm3, species %in% transients3) %>% 
  group_by(period) %>% 
  summarise(transients = sum(abundance))
rel_abund3 <- inner_join(rel_abund3, total_by_period3)
rel_abund_transients3 <- mutate(rel_abund3, transient_rel_abund = transients/total)

# add date to abundances
for (i in 1:length(rdat$month)){
  if (str_length(rdat$month[i]) == 1){
    rdat$month[i] = paste('0', rdat$month[i], sep = '')
  }
}

dates <- dplyr::select(rdat, year, month, period) %>% 
  tidyr::unite(date, year, month, sep = '-') %>% 
  distinct()
dates$period <- as.integer(dates$period)

abund_dates1 <- inner_join(rel_abund_transients1, dates, by = 'period') 
abund_dates2 <- inner_join(rel_abund_transients2, dates, by = 'period') 
abund_dates3 <- inner_join(rel_abund_transients3, dates, by = 'period') 
abund_dates <- bind_rows(abund_dates1, abund_dates2, abund_dates3)

# ADD PLOTS TRAPPED #

# get data on how many plots were trapped each period
plots_trapped <- trapping %>% 
  dplyr::select(period, plot, sampled) %>% 
  group_by(period) %>% 
  summarise(count_plots = sum(sampled))

abund_dates <- left_join(abund_dates, plots_trapped, by = "period") 
abund_dates <- filter(abund_dates, date < 2015)

# if number of plots trapped is between 10 and 17, double the transients
for (i in 1:nrow(abund_dates)){
  if (abund_dates$count_plots[i] <= 17){
    abund_dates$transient_rel_abund[i] <- abund_dates$transient_rel_abund[i]*2
  }
}


### NDVI ###

# make a date column for NDVI
NDVI <- arrange(NDVI, year, month) %>% 
  rename(ndvi = x)

for (i in 1:nrow(NDVI)){
  if (str_length(NDVI$month[i]) == 1){
    NDVI$month[i] = paste('0', NDVI$month[i], sep = '')
  }
}

NDVI <- dplyr::select(NDVI, year, month, ndvi) %>% 
  tidyr::unite(date, year, month, sep = '-', remove = FALSE) %>% 
  distinct()

# make a column for NDVI centered around the median
NDVI_peak <- mutate(NDVI, NDVIpeak = ndvi - median(ndvi, na.rm = TRUE))

# filter NDVI by community dates
NDVI1 <- filter(NDVI_peak, date %in% unique(abund_dates1$date))
NDVI2 <- filter(NDVI_peak, date %in% unique(abund_dates2$date))
NDVI3 <- filter(NDVI_peak, date %in% unique(abund_dates3$date))
NDVI_all <- bind_rows(NDVI1, NDVI2, NDVI3)
NDVI_all <- filter(NDVI_all, date < 2015)
#         maybe include December of 1988 and January/February of 2015
#         so the seasonal sums/means are more accurate?

### COMBINE TRANSIENTS AND NDVI ###

combined_data <- full_join(abund_dates, NDVI_all, by = "date")
combined_data$date = as.yearmon(combined_data$date)

# deal with periods that have two different months
# select first month out of the two?

combined_data <- combined_data[!duplicated(combined_data$period),]

# get all yearmon dates
all_dates <- as.data.frame(seq(as.yearmon(min(NDVI_peak$date)), as.yearmon(max(NDVI_all$date)), by = 1/12))
colnames(all_dates) = c("date")
#additional_year <- as.data.frame(seq(from = max(all_dates$date) + 1/12, by = 1/12, length.out = 12))
#colnames(additional_year) = c("date")
#all_dates <- rbind(all_dates, additional_year)

combined_data <- full_join(combined_data, all_dates, by = "date") %>% 
  arrange(date)
combined_data <- combined_data[-c(1:74),]
combined_data$month <- as.integer(combined_data$month)

# add in NDVI values even if there aren't transient values
for (i in 1:nrow(combined_data)){
  
  # Is the value NA? If so, save the date as object 'date'
  if (is.na(combined_data$ndvi[i])){
    date = combined_data$date[i]
    
    # Is it true that the date is in the NDVI_all dataframe?
    # If so, replace the NAs in the ndvi column and NDVIpeak columns with values
    #     from the NDVI_all dataframe
    if (length(which(NDVI_all$date == date))>0){
      combined_data$ndvi[i] <- NDVI_all[NDVI_all$date == date, "ndvi"]
      combined_data$NDVIpeak[i] <- NDVI_all[NDVI_all$date == date, "NDVIpeak"]
      
    }
  }
}

# WITHOUT INTERPOLATION # 

# fill in any missing months and years
for (i in 2:nrow(combined_data)){
  
  # fill in missing month and year values
  if (is.na(combined_data$month[i])){
    combined_data$year[i] <- year(combined_data$date[i])
    combined_data$month[i] <- month(combined_data$date[i])
  }
}

# combined_data$year[138:139] <- 1999
# combined_data$month[138] <- 7
# combined_data$month[139] <- 8


# ADD SEASONS #

seasonal_data <- portalr::add_seasons(combined_data, season_level = 4)

# add seasonal year column for summary
seasonal_data$seasonal_year <- NA
for (i in 1:nrow(seasonal_data)){
  if(month(seasonal_data$date[i]) == 12){
    seasonal_data$seasonal_year[i] <- (year(seasonal_data$date[i]) + 1)
  } else {
    seasonal_data$seasonal_year[i] <- year(seasonal_data$date[i])
  }
}

# summarize by season per year
# seasonal_data <- seasonal_data[-c(1:40),]
by_season <- seasonal_data %>% 
  group_by(seasonal_year, season) %>% 
  summarise(transients = mean(transient_rel_abund, na.rm = TRUE), ndvi = mean(NDVIpeak, na.rm = TRUE)) 
by_season <- by_season[-c(1:5),]

# make into quarterly data for easier plotting
by_season$quarter <- NA
for (i in 1:nrow(by_season)){
  if(by_season$season[i] == 'winter'){
    by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 1)
  } else if (by_season$season[i] == 'spring'){
    by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 2)
  } else if (by_season$season[i] == 'summer'){
    by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 3)
  } else {
    by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 4)
  }  
  
}

by_season$quarter <- as.yearqtr(by_season$quarter)

# plot the time series
ggplot(by_season) +
  geom_point(aes(x = quarter, y = ndvi), color = "dark green") +
  geom_line(aes(x = quarter, y = ndvi), color = "dark green") +
  geom_point(aes(x = quarter, y = transients), color = "salmon") +
  geom_line(aes(x = quarter, y = transients), color = "salmon")

# run correlation functions

by_season_NAomit <- arrange(by_season, quarter) %>% 
  na.omit()

acf(by_season_NAomit$ndvi)  
acf(by_season_NAomit$transients)
ccf(by_season_NAomit$transients, by_season$ndvi)


# WITH INTERPOLATION #

for (i in 2:nrow(combined_data)){
  
  # fill in missing transient values
  if (is.na(combined_data$transient_rel_abund[i])){
    combined_data$transient_rel_abund[i] <- (combined_data$transient_rel_abund[i-1] + combined_data$transient_rel_abund[i+1])/2
  }
  
  # fill in missing NDVI values
  if (is.na(combined_data$NDVIpeak[i])){
    if(!is.na(combined_data$NDVIpeak[i-1]))
    combined_data$NDVIpeak[i] <- (combined_data$NDVIpeak[i-1] + combined_data$NDVIpeak[i+1])/2
  } 
  
  # fill in missing month and year values
  if (is.na(combined_data$month[i])){
    combined_data$year[i] <- year(combined_data$date[i])
    combined_data$month[i] <- month(combined_data$date[i])
  }
  
}

combined_data$transient_rel_abund[212:213] <- 0
combined_data$NDVIpeak[212:213] <- (combined_data$NDVIpeak[211] + combined_data$NDVIpeak[214])/2

# ADD SEASONS #

seasonal_data <- portalr::add_seasons(combined_data, season_level = 4)

# add seasonal year column for summary
seasonal_data$seasonal_year <- NA
for (i in 1:nrow(seasonal_data)){
  if(seasonal_data$month[i] == 12){
    seasonal_data$seasonal_year[i] <- (seasonal_data$year[i] + 1)
  } else {
    seasonal_data$seasonal_year[i] <- seasonal_data$year[i]
  }
}

# summarize by season per year
seasonal_data <- seasonal_data[-c(1:40),]
by_season <- seasonal_data %>% 
  group_by(seasonal_year, season) %>% 
  summarise(transients = mean(transient_rel_abund, na.rm = TRUE), ndvi = mean(NDVIpeak, na.rm = TRUE)) 

# make into quarterly data for easier plotting
by_season$quarter <- NA
for (i in 1:nrow(by_season)){
  if(by_season$season[i] == 'winter'){
    by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 1)
  } else if (by_season$season[i] == 'spring'){
      by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 2)
  } else if (by_season$season[i] == 'summer'){
    by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 3)
  } else {
      by_season$quarter[i] = paste0(by_season$seasonal_year[i], "-", 4)
  }  
  
}
  
by_season$quarter <- as.yearqtr(by_season$quarter)
by_season <- by_season[-92,]

# plot the time series
ggplot(by_season) +
  geom_point(aes(x = quarter, y = ndvi), color = "dark green") +
  geom_line(aes(x = quarter, y = ndvi), color = "dark green") +
  geom_point(aes(x = quarter, y = transients), color = "salmon") +
  geom_line(aes(x = quarter, y = transients), color = "salmon")

# try running acfs and ccfs

by_season <- arrange(by_season, quarter)

acf(by_season$ndvi)  
acf(by_season$transients)
ccf(by_season$transients, by_season$ndvi)
