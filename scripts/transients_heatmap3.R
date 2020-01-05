# Transients w/ Portalr (Lunar)
  # most code taken from transients_heatmap2 and adjusted for newmoons
# EKB
# Jan 2020

#=========================================================== 
# LIBRARIES and DATA #

library(RCurl)
library(tidyverse)
library(portalr)

trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
moon_dates <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"))


#=========================================================== 
# FUNCTIONS #

find_reg_persist <- function(data){
  # create column for relative regional persistence (periods overall)
  
  # get summed abundance by species and period
  dat <- data %>% 
    dplyr::group_by(species, newmoonnumber) %>% 
    dplyr::summarise(abundance = sum(abundance, na.rm = TRUE))
  
  # get summed number of periods each species is present
  dat1 <- dplyr::select(dat, species) %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(newmoonnumber = n())
  
  total_periods <- length(unique(data$newmoonnumber))
  
  # calculate regional persistance column
  dat1 <- dplyr::mutate(dat1, rel_reg_persist = newmoonnumber/total_periods)
  
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
    group_by(newmoonnumber) %>% 
    summarise(total = sum(abundance))
  
  # get relative abundance of all transients per period
  rel_abund_trans <- data %>% 
    filter(species %in% transients) %>% 
    group_by(newmoonnumber) %>% 
    summarise(transients = sum(abundance)) %>% 
    inner_join(total_by_period) %>% 
    mutate(trans_rel_abund = transients/total)
  
  return(rel_abund_trans)
  
}
