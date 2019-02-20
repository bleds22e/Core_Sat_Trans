# Transients Heatmaps 2
# EKB
# February 2019

# LIBRARIES and DATA #

library(RCurl)
library(zoo)
library(tidyverse)

rdat <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"), na.strings = '')
trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))
  # different NDVI file

# FUNCTIONS #

find_reg_persist <- function(data){
  # create column for relative regional persistence (periods overall)
  
  # get summed abundance by species and period
  dat <- data %>% 
    dplyr::group_by(species, period) %>% 
    dplyr::summarise(abundance = sum(abundance, na.rm = TRUE))
  
  # get summed number of periods each species is present
  dat1 <- dplyr::select(dat, species) %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(periods = n())
  
  total_periods <- length(unique(data$period))
  
  # calculate regional persistance column
  dat1 <- dplyr::mutate(dat1, rel_reg_persist = periods/total_periods)
  
  return(dat1)
  
}


get_transient_rel_abund <- function(data){
  # select transient species and get total transient relative abundance
  
  # get list of transient species
  transients <- data %>% 
    filter(abundance > 0) %>% 
    find_reg_persist() %>% 
    filter(rel_reg_persist <= 1/3, species != 'PB') %>% 
    select(species)
  transients <- unlist(apply(transients, 1, list), recursive = FALSE)
  
  # get total abundance by period
  total_by_period <- data %>% 
    group_by(period) %>% 
    summarise(total = sum(abundance))
  
  # get relative abundance of all transients per period
  rel_abund_trans <- data %>% 
    filter(species %in% transients) %>% 
    group_by(period) %>% 
    summarise(transients = sum(abundance)) %>% 
    inner_join(total_by_period) %>% 
    mutate(trans_rel_abund = transients/total)
  
  return(rel_abund_trans)
  
}

#===============================================================
# GET TRANSIENT RODENT ABUNDANCES 


# get rodent abundance data
abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                            plots = "all", unknowns = F, shape = "flat", 
                            time = c("period"), min_plots = 10)

# pre 1984, 1984-2000, 2000-2009, post-2009 

comm1 <- filter(abund, period <= 75) # pre-1984
comm2 <- filter(abund, period > 75 & period < 261) # 1984-2000
comm3 <- filter(abund, period >= 261 & period < 381) # 2000-2009
comm4 <- filter(abund, period >= 381) # post-2009

# get transient relative abundance
rel_abund_trans1 <- get_transient_rel_abund(comm1)
rel_abund_trans2 <- get_transient_rel_abund(comm2)
rel_abund_trans3 <- get_transient_rel_abund(comm3)
rel_abund_trans4 <- get_transient_rel_abund(comm4)

rel_abund_trans <- do.call("rbind", list(rel_abund_trans1,
                                         rel_abund_trans2,
                                         rel_abund_trans3,
                                         rel_abund_trans4))

# get year and month columns paired with period column
trapping_dates <- trapping %>% 
  select(month, year, period) %>% 
  group_by(period) %>% 
  unique()

# get number of plots trapped per period
trapping_effort <- trapping %>% 
  select(period, sampled) %>% 
  group_by(period) %>% 
  summarise(nplots = sum(sampled)) 

# add month, year, and effor to relative transient abundance
trapping_data <- full_join(trapping_dates, trapping_effort)
rel_abund_data <- full_join(trapping_data, rel_abund_trans)

# double the number of transients if approriate
for (i in 1:nrow(rel_abund_data)){
  if (rel_abund_data$nplots[i] <= 17){
    rel_abund_data$transient_rel_abund[i] <- rel_abund_data$transient_rel_abund[i]*2
  }
}

#=========================================================
# COMBINE NDVI AND TRANSIENT DATA

# put months in order
NDVI <- NDVI %>% arrange(year, month)

# add transient relative abundance and yearmon object
NDVI_transient <- full_join(rel_abund_data, NDVI)
NDVI_transient$date <- as.yearmon(paste(NDVI_transient$year, NDVI_transient$month), "%Y %m")  

################### need to deal with multiple months per period 
                  # and multiple periods per month
                  # would using new_moons as the x-axis fix this issue?

#==========================================================
# MAKE NEW DATAFRAME FOR HEATMAP #

# get median value of NDVI
NDVI_transient <- mutate(NDVI_transient, NDVI_median_center = ndvi - median(ndvi, na.rm = TRUE))

# make outline of data frame
nrows <- length(NDVI_transient$date)
lags <- data.frame(NDVI = NDVI_transient$NDVI_median_center,
                   d.NDVI = numeric(nrows),
                   transientT0 = numeric(nrows),
                   transientT1 = numeric(nrows),
                   transientT2 = numeric(nrows),
                   transientT3 = numeric(nrows),
                   transientT4 = numeric(nrows),
                   transientT5 = numeric(nrows),
                   transientT6 = numeric(nrows),
                   transientT7 = numeric(nrows),
                   transientT8 = numeric(nrows),
                   transientT9 = numeric(nrows),
                   transientT10 = numeric(nrows),
                   transientT11 = numeric(nrows))

for (i in 2:length(combined_data$date)){
  
  # find change in NDVI
  lags$d.NDVI[i] <- lags$NDVI[i] - lags$NDVI[i-1]
  
  # find the transient rel abundances
  lags$transientT0[i] <- combined_data$transient_rel_abund[i]
  lags$transientT1[i] <- combined_data$transient_rel_abund[i + 1]
  lags$transientT2[i] <- combined_data$transient_rel_abund[i + 2]
  lags$transientT3[i] <- combined_data$transient_rel_abund[i + 3]
  lags$transientT4[i] <- combined_data$transient_rel_abund[i + 4]
  lags$transientT5[i] <- combined_data$transient_rel_abund[i + 5]
  lags$transientT6[i] <- combined_data$transient_rel_abund[i + 6]
  lags$transientT7[i] <- combined_data$transient_rel_abund[i + 7]
  lags$transientT8[i] <- combined_data$transient_rel_abund[i + 8]
  lags$transientT9[i] <- combined_data$transient_rel_abund[i + 9]
  lags$transientT10[i] <- combined_data$transient_rel_abund[i + 10]
  lags$transientT11[i] <- combined_data$transient_rel_abund[i + 11]
  
}