# Transients & NDVI by species
# EKB
# March 2019

# LIBRARIES and DATA #

library(RCurl)
library(zoo)
library(tidyverse)

rdat <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"), na.strings = '')
trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))

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


get_transient_sp_rel_abund <- function(data){
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
    group_by(period, species) %>% 
    inner_join(total_by_period) %>% 
    mutate(trans_rel_abund = abundance/total)
  
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
rel_abund_trans_sp1 <- get_transient_sp_rel_abund(comm1)
rel_abund_trans_sp2 <- get_transient_sp_rel_abund(comm2)
rel_abund_trans_sp3 <- get_transient_sp_rel_abund(comm3)
rel_abund_trans_sp4 <- get_transient_sp_rel_abund(comm4)

rel_abund_trans_sp <- do.call("rbind", list(rel_abund_trans_sp1,
                                         rel_abund_trans_sp2,
                                         rel_abund_trans_sp3,
                                         rel_abund_trans_sp4))

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
rel_abund_data <- full_join(trapping_data, rel_abund_trans_sp)
