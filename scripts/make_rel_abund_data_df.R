# script to make rel_abund_data dataframe for making adjustments to transients 
# by species dataframes

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
    rel_abund_data$trans_rel_abund[i] <- rel_abund_data$trans_rel_abund[i]*2
  }
}

### nudge some periods around then select second period w/ same date

rel_abund_data <- rel_abund_data %>% 
  arrange(year, month)

# select certain months of periods with 2 months
rel_abund_data[83, 1] <- 6
rel_abund_data[112:113, 3:7] <- rel_abund_data[113:114, 3:7]
rel_abund_data[143, 1] <- 8
rel_abund_data[178, 1] <- 6
rel_abund_data[197, 1] <- 1
rel_abund_data[214, 1] <- 5
rel_abund_data[220, 1] <- 11
rel_abund_data[240, 1] <- 8
rel_abund_data[247, 1] <- 2
rel_abund_data[249, 1] <- 4
rel_abund_data[272, 1] <- 5
rel_abund_data[273, 1] <- 6
rel_abund_data[281, 1] <- 2
rel_abund_data[303, 1:2] <- c(1, 2003)
rel_abund_data[307, 1] <- 4
rel_abund_data[342, 1] <- 2
rel_abund_data[403, 1] <- 6
rel_abund_data[408, 1] <- 11

rel_abund_data <- rel_abund_data[-c(48, 52, 56, 85, 114, 118, 148, 179, 212, 
                                    246, 276, 305, 336, 340, 370, 405),]