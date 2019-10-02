# Transients & NDVI by species
# EKB
# March 2019

# LIBRARIES and DATA #

library(RCurl)
library(zoo)
library(tidyverse)

trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))

source("scripts/make_rel_abund_data_df.R")

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
comm2 <- filter(abund, period > 75 & period < 261) # 1984-1999
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
trapping_dates_sp <- trapping %>% 
  select(month, year, period) %>% 
  group_by(period) %>% 
  unique()

# get number of plots trapped per period
trapping_effort_sp <- trapping %>% 
  select(period, sampled) %>% 
  group_by(period) %>% 
  summarise(nplots = sum(sampled)) 

# add month, year, and effot to relative transient abundance
trapping_data_sp <- full_join(trapping_dates_sp, trapping_effort_sp)
rel_abund_data_sp <- full_join(trapping_data_sp, rel_abund_trans_sp) %>% 
  filter(nplots >= 20) # remove periods with fewer than 20 plots

# add in NDVI data and nudge periods
NDVI_transient <- full_join(rel_abund_data_sp, NDVI)

NDVI_transient <- NDVI_transient %>% 
  ungroup() %>% 
  select(-period)
rel_abund_data <- rel_abund_data %>% 
  select(-(nplots:trans_rel_abund))
NDVI_transient <- full_join(rel_abund_data, NDVI_transient)

# make sure there's a row for every yearmon combo
all_dates <- expand.grid(month = seq(1:12), year = 1978:2018)
NDVI_transient <- full_join(NDVI_transient, all_dates, by = c("month" = "month", "year" = "year"))

# put months in order
NDVI_transient <- NDVI_transient %>% arrange(year, month)

# fill in easily interpolated missing values
for (i in 2:nrow(NDVI_transient)){
  
  # fill in missing transient values
  if (is.na(NDVI_transient$trans_rel_abund[i])){
    NDVI_transient$trans_rel_abund[i] <- (NDVI_transient$trans_rel_abund[i-1] + NDVI_transient$trans_rel_abund[i+1])/2
  }
  
  # fill in missing NDVI values
  if (is.na(NDVI_transient$ndvi[i])){
    if(!is.na(NDVI_transient$ndvi[i-1]))
      NDVI_transient$ndvi[i] <- (NDVI_transient$ndvi[i-1] + NDVI_transient$ndvi[i+1])/2
  } 
  
}

# add in rows for each species with NAs based on community

NDVI_transient_NA <- NDVI_transient[is.na(NDVI_transient$species),]

trans1 <- as.data.frame(unique(rel_abund_trans_sp1$species)) 
colnames(trans1) <- c("species")
trans2 <- as.data.frame(unique(rel_abund_trans_sp2$species)) 
colnames(trans2) <- c("species")
trans3 <- as.data.frame(unique(rel_abund_trans_sp3$species)) 
colnames(trans3) <- c("species")
trans4 <- as.data.frame(unique(rel_abund_trans_sp4$species)) 
colnames(trans4) <- c("species")

comm1_na <- NDVI_transient_NA[NDVI_transient_NA$year < 1984,]
comm2_na <- filter(NDVI_transient_NA, year > 1983 & year < 2000)
comm3_na <- filter(NDVI_transient_NA, year >= 2000 & year < 2010)
comm4_na <- filter(NDVI_transient_NA, year >= 2010)

comm_na_species <- list()
n = 1

communities <- list(comm1_na, comm2_na, comm3_na, comm4_na)
transient_sp <- list(trans1, trans2, trans3, trans4)

for(c in 1:length(communities)){
  
  comm <- communities[[c]]
  trans <- transient_sp[[c]]
  
  for (i in 1:nrow(comm)){
    
    row <- comm[i,]
    
    for(p in 1:nrow(trans)){
      
      row$species = trans[p,]
      comm_na_species[[n]] <- row
      n = n +1
      
    }
  }
}

comm_na_species <- do.call(rbind, comm_na_species)

for(i in 1:nrow(NDVI_transient)){
  if(is.na(NDVI_transient$period[i]) & !is.na(NDVI_transient$species[i])){
    NDVI_transient$trans_rel_abund[i] <- NA
  }
}


NDVI_transient <- NDVI_transient[!is.na(NDVI_transient$species),]
NDVI_transient <- rbind(NDVI_transient, comm_na_species)

#==========================================================
# PLOTTING 
# plot NDVI vs transients by species -- time series

# make ndvi median value and make long df
NDVI_transient$ndvi <- NDVI_transient$ndvi - median(NDVI_transient$ndvi, na.rm = TRUE)
NDVI_transient_long <- gather(NDVI_transient, key = "data_type", value = "value", 8:9)

# remove rows in species = NA
NDVI_transient_long <- NDVI_transient_long[!is.na(NDVI_transient_long$species),]

# Time Series of NDVI and transients by species #
ggplot(data = NDVI_transient_long, aes(x = period, y = value)) +
  geom_line(aes(color = data_type), size = 1) +
  facet_wrap(~ species, ncol = 3)
ggsave("plots/GIMMs_plots/timeseries_by_species.png")

#==========================================================
# MAKE NEW DATAFRAME FOR HEATMAP #

NDVI_transient$date <- as.yearmon(paste(NDVI_transient$year, NDVI_transient$month), "%Y %m") 
species_dfs_list <- split(NDVI_transient, NDVI_transient$species)
species_dfs_list$PH1 <- filter(species_dfs_list$PH, year < 2000)
species_dfs_list$PH2 <- filter(species_dfs_list$PH, year > 2009)
species_dfs_list$SH1 <- filter(species_dfs_list$SH, year < 2000)
species_dfs_list$SH2 <- filter(species_dfs_list$SH, year > 2009)
species_dfs_list$PH <- NULL
species_dfs_list$SH <- NULL

species_dfs_list <- species_dfs_list[sapply(species_dfs_list, function(x) dim(x)[1] > 0)]

# make outline of data frame
lags_dfs <- list()

for (df in 1:length(species_dfs_list)){
  
  nrows <- nrow(species_dfs_list[[df]])
  lags <- data.frame(species = character(nrows),
                     NDVI = numeric(nrows),
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
  
  lags$species <- species_dfs_list[[df]]$species
  lags$NDVI <- species_dfs_list[[df]]$ndvi
  
  for (i in 2:nrows){
    
    # find change in NDVI
    lags$d.NDVI[i] <- lags$NDVI[i] - lags$NDVI[i-1]
    
    # find the transient rel abundances
    lags$transientT0[i] <- species_dfs_list[[df]]$trans_rel_abund[i]
    lags$transientT1[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 1]
    lags$transientT2[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 2]
    lags$transientT3[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 3]
    lags$transientT4[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 4]
    lags$transientT5[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 5]
    lags$transientT6[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 6]
    lags$transientT7[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 7]
    lags$transientT8[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 8]
    lags$transientT9[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 9]
    lags$transientT10[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 10]
    lags$transientT11[i] <- species_dfs_list[[df]]$trans_rel_abund[i + 11]
    
  }
  
  lags_dfs[[df]] <- lags
  
}

lags_dfs <- do.call(rbind, lags_dfs)
lags_dfs <- lags_dfs[complete.cases(lags_dfs),]

lags_long <- tidyr::gather(lags_dfs, "time_lag", "transients", 4:15)

new_levels = c("transientT0", "transientT1", "transientT2",
               "transientT3", "transientT4", "transientT5",
               "transientT6", "transientT7", "transientT8",
               "transientT9", "transientT10", "transientT11")
lags_long <- arrange(mutate(lags_long, time_lag = factor(time_lag, levels = new_levels)), time_lag)


### Get "contingency" tables ###

lags_dfs_list <- split(lags_long, lags_long$species)
lags_dfs_list <- lags_dfs_list[sapply(lags_dfs_list, function(x) dim(x)[1] > 0)]

quadrant_means_list <- list()

for (sp in 1:length(lags_dfs_list)){
  
  data <- lags_dfs_list[[sp]]
  
  x.neg_y.pos <- data %>% 
    filter(NDVI < 0, d.NDVI > 0, transients > 0.005) %>% 
    group_by(time_lag) %>% 
    summarise(x.neg_y.pos = round(mean(transients), 5))
  x.pos_y.pos <- data %>% 
    filter(NDVI >= 0, d.NDVI > 0, transients > 0.005) %>% 
    group_by(time_lag) %>% 
    summarise(x.pos_y.pos = round(mean(transients), 5))
  x.pos_y.neg <- data %>% 
    filter(NDVI >= 0, d.NDVI <= 0, transients > 0.005) %>% 
    group_by(time_lag) %>% 
    summarise(x.pos_y.neg = round(mean(transients), 5))
  x.neg_y.neg <- data %>% 
    filter(NDVI < 0, d.NDVI <= 0, transients > 0.005) %>% 
    group_by(time_lag) %>% 
    summarise(x.neg_y.neg = round(mean(transients), 5))
  
  quadrant_means <- plyr::join_all(list(x.neg_y.pos, x.pos_y.pos, x.pos_y.neg, x.neg_y.neg), by = "time_lag")
  quadrant_means_long <- gather(quadrant_means, key = "quadrant", value = "mean_transients", 2:5)
  quadrant_means_long$species <- lags_dfs_list[[sp]]$species[1]
  quadrant_means_list[[sp]] <- quadrant_means_long
  
}

quadrant_means_by_sp <- do.call(rbind, quadrant_means_list)
quadrant_means_by_sp <- quadrant_means_by_sp[complete.cases(quadrant_means_by_sp),]

ggplot(quadrant_means_by_sp, aes(time_lag, mean_transients, color = quadrant, group = quadrant)) + 
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ species) +
  theme_bw() +
  xlab("Lag Time") +
  ylab("Mean Transients (above 0.005)") +
  theme(axis.text.x = element_text(angle = -45, hjust = -.1))
ggsave("plots/GIMMs_plots/species_through_time_by_quadrant_0.005.png")


#============================================================
# WORK SPACE


