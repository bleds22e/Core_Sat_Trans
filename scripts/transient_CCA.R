# Transient CCA
# EKB
# Oct 2018

# LIBRARIES #
library(RCurl)
library(tidyverse)
library(stringr)
library(zoo)

# FUNCTION #

find_reg_persist <- function(data){
  # create column for relative regional persistence (periods overall)
  dat <- select(data, species, period) %>% 
    group_by(species, period) %>% 
    summarise(period_count = n())
  dat1 <- select(dat, species) %>% 
    group_by(species) %>% 
    summarise(periods = n())
  total_periods <- length(unique(data$period))
  dat1 <- mutate(dat1, rel_reg_persist = periods/total_periods)
  return(dat1)
}

# GET DATA #

rdat <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"), na.strings = '')
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))

#===============================================================

abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                            plots = "all", unknowns = F, shape = "flat", time = "period")

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

# par(mfrow = c(3,1))
# hist(comm1_persist$rel_reg_persist)
# hist(comm2_persist$rel_reg_persist)
# hist(comm3_persist$rel_reg_persist)

#dev.copy(png, "plots/hist_portal_comms.png")
#dev.off()

# transients for each community

transients1 <- comm1_persist %>% 
  filter(rel_reg_persist < 0.5, species != "PB") %>% 
  select(species)
transients1 <- unlist(apply(transients1, 1, list), recursive = FALSE)

transients2 <- comm2_persist %>% 
  filter(rel_reg_persist < 0.5) %>% 
  select(species)
transients2 <- unlist(apply(transients2, 1, list), recursive = FALSE)

transients3 <- comm3_persist %>% 
  filter(rel_reg_persist < 0.5) %>% 
  select(species)
transients3 <- unlist(apply(transients3, 1, list), recursive = FALSE)


# GET ABUNDANCES #

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

dates <- select(rdat, year, month, period) %>% 
  tidyr::unite(date, year, month, sep = '-') %>% 
  distinct()
dates$period <- as.integer(dates$period)

# make a date column for NDVI
NDVI <- arrange(NDVI, year, month)

for (i in 1:length(NDVI$month)){
  if (str_length(NDVI$month[i]) == 1){
    NDVI$month[i] = paste('0', NDVI$month[i], sep = '')
  }
}

NDVI <- select(NDVI, year, month, ndvi) %>% 
  tidyr::unite(date, year, month, sep = '-') %>% 
  distinct()

# filter NDVI by community dates
abund_dates1 <- inner_join(rel_abund_transients1, dates, by = 'period') 
abund_dates2 <- inner_join(rel_abund_transients2, dates, by = 'period') 
abund_dates3 <- inner_join(rel_abund_transients3, dates, by = 'period') 
abund_dates <- bind_rows(abund_dates1, abund_dates2, abund_dates3)

NDVI_peak <- mutate(NDVI, NDVIpeak = ndvi - median(ndvi))

NDVI1 <- filter(NDVI_peak, date %in% unique(abund_dates1$date))
NDVI2 <- filter(NDVI_peak, date %in% unique(abund_dates2$date))
NDVI3 <- filter(NDVI_peak, date %in% unique(abund_dates3$date))
NDVI_all <- bind_rows(NDVI1, NDVI2, NDVI3)
NDVI_all$date = as.yearmon(NDVI_all$date)

combined_data <- inner_join(abund_dates, NDVI_all, by = "date")
combined_data$date = as.yearmon(combined_data$date)

# deal with periods that have two different months
# select first month out of the two?

combined_data <- combined_data[!duplicated(combined_data$period),]

# get all yearmon dates
all_dates <- as.data.frame(seq(as.yearmon(min(NDVI_peak$date)), as.yearmon(max(NDVI_peak$date)), by = 1/12))
colnames(all_dates) = c("date")

combined_data <- full_join(all_dates, combined_data, by = "date")
combined_data <- combined_data[-c(1:30),]

# quantiles
upper <- as.numeric(quantile(NDVI_peak$NDVIpeak, .75))
lower <- as.numeric(quantile(NDVI_peak$NDVIpeak, .25))

# add in NDVI values even if there aren't transient values
for (i in 1:length(combined_data$date)){
  
  # Is the value NA? If so, save the date as object 'date'
  if (is.na(combined_data$ndvi[i])){
    date = combined_data$date[i]
    
    # Is it true that the date is in the NDVI_all dataframe?
    # If so, replace the NAs in the ndvi column and NDVIpeak columns with values
    # from the NDVI_all dataframe
    if (length(which(NDVI_all$date == date))>0){
      combined_data$ndvi[i] <- NDVI_all[NDVI_all$date == date, "ndvi"]
      combined_data$NDVIpeak[i] <- NDVI_all[NDVI_all$date == date, "NDVIpeak"]
      
    }
  }
}

combined_data$NDVI_quantiles = NA

for (i in 1:length(combined_data$NDVIpeak)){
  if (!is.na(combined_data$NDVIpeak[i])){
    if (combined_data$NDVIpeak[i] > upper){
      combined_data$NDVI_quantiles[i] = combined_data$NDVIpeak[i]
    } else {
      combined_data$NDVI_quantiles[i] = NA
    }
  }
}




# first try to get pulse ID
# then can work from there fairly easily, I think
combined_data$pulseID <- NA
pulseID = 1

for (i in 1:length(combined_data$date)){
  
  if (!is.na(combined_data$NDVI_quantiles[i])){
    if(is.na(combined_data$pulseID[i-1])){
      combined_data$pulseID[i] <- pulseID
      pulseID = pulseID + 1
    } else{
      combined_data$pulseID[i] <- combined_data$pulseID[i-1]
    }
    
  }
  
}


# make outline of data frame
pulses <- data.frame(pulseID = integer(),
                     pulseMax = numeric(),
                     pulseDuration = integer(),
                     pulseSum = numeric(),
                     transientT1 = numeric(),
                     tranisentT2 = numeric(),
                     transientT3 = numeric(),
                     transientT4 = numeric())


#====================================================================
# Working Area #

