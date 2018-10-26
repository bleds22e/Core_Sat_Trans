# Transient CCA
# EKB
# Oct 2018

# LIBRARIES #
library(RCurl)
library(stringr)
library(zoo)
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
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))

#===============================================================

abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                            plots = "all", unknowns = F, shape = "flat", time = "period")

### ALERT ###

# missing some periods because not all periods were fully trapped
# also missing some NDVI values because of...reasons
# makes for *lots* of NA values, making CCA analysis impossible
# should I:
#   - for periods where site was only trapped for one night,
#     double the number of transients caught while trapping?
#   - for NDVI, take the average of before and after missing value?
#   - what if there are 2 NDVI missing in a row?




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

combined_data <- full_join(abund_dates, NDVI_all, by = "date")
combined_data$date = as.yearmon(combined_data$date)

# deal with periods that have two different months
# select first month out of the two?

combined_data <- combined_data[!duplicated(combined_data$period),]

# get all yearmon dates
all_dates <- as.data.frame(seq(as.yearmon(min(NDVI_peak$date)), as.yearmon(max(NDVI_peak$date)), by = 1/12))
colnames(all_dates) = c("date")

combined_data <- full_join(combined_data, all_dates, by = "date")
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

# get pulseID
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

# make data frame of rows only with pulses
combined_pulses <- combined_data[!is.na(combined_data$pulseID),]

# make outline of data frame
nrows <- max(combined_pulses$pulseID, na.rm = TRUE)
pulses <- data.frame(pulseID = integer(nrows),
                     pulseMax = numeric(nrows),
                     pulseDuration = integer(nrows),
                     pulseSum = numeric(nrows),
                     timeBetween = integer(nrows), # time from end of peak to start of next
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
                     transientT11 = numeric(nrows),
                     transientT12 = numeric(nrows))

pulse_seq <- seq(from = 1, to = nrows)
pulses$pulseID <- pulse_seq

for (i in pulse_seq){
  
  # select only rows from that pulse
  current_pulse <- combined_pulses[combined_pulses$pulseID == i,]
  
  # find the maximum NDVI value
  pulses$pulseMax[i] <- max(current_pulse$NDVIpeak)
  
  # find pulse duration
  pulses$pulseDuration[i] <- nrow(current_pulse)
  
  # find pulse sum
  pulses$pulseSum[i] <- sum(current_pulse$NDVIpeak)
  
  # find time (months) between end of current pulse and start of next pulse
  pulses$timeBetween[i] <- 12 * (min(combined_pulses[combined_pulses$pulseID == i + 1, "date"]) - max(current_pulse$date))
  
  # find the transient rel abundances
  min_date <- min(current_pulse$date)
  min_row <- as.integer(rownames(combined_data[combined_data$date == min_date,]))
  pulses$transientT0[i] <- combined_data[as.character(min_row), "transient_rel_abund"]
  pulses$transientT1[i] <- combined_data[as.character(min_row + 1), "transient_rel_abund"]
  pulses$transientT2[i] <- combined_data[as.character(min_row + 2), "transient_rel_abund"]
  pulses$transientT3[i] <- combined_data[as.character(min_row + 3), "transient_rel_abund"]
  pulses$transientT4[i] <- combined_data[as.character(min_row + 4), "transient_rel_abund"]
  pulses$transientT5[i] <- combined_data[as.character(min_row + 5), "transient_rel_abund"]
  pulses$transientT6[i] <- combined_data[as.character(min_row + 6), "transient_rel_abund"]
  pulses$transientT7[i] <- combined_data[as.character(min_row + 7), "transient_rel_abund"]
  pulses$transientT8[i] <- combined_data[as.character(min_row + 8), "transient_rel_abund"]
  pulses$transientT9[i] <- combined_data[as.character(min_row + 9), "transient_rel_abund"]
  pulses$transientT10[i] <- combined_data[as.character(min_row + 10), "transient_rel_abund"]
  pulses$transientT11[i] <- combined_data[as.character(min_row + 11), "transient_rel_abund"]
  pulses$transientT12[i] <- combined_data[as.character(min_row + 12), "transient_rel_abund"]
  
}

pulses$timeBetween[38] <- NA

#========================================================================
# Run CCA Analysis

# separate into env and response variables
pulses_NDVI <- pulses[,1:5]
pulses_transients <- pulses[,6:18]

# histograms of transient distributions
mapply(hist, as.data.frame(pulses_transients[,1:13], 
                           main = colnames(pulses_transients[,1:13]),
                           xlab = "abundance"))

# log transformation
log.full <- log1p(pulses_transients)

# check row and column sum variability
rsum <- rowSums(log.full, na.rm = TRUE)
csum <- colSums(log.full, na.rm = TRUE)
hist(rsum)
hist(csum)
cv(rsum)
cv(csum)

# standardize the rows because cv > 50
rTrans <- sweep(log.full, 1, rsum, "/")
rTrans <- rTrans[-c(34,35),]
cv(rowSums(rTrans, na.rm = TRUE))
cv(colSums(rTrans, na.rm = TRUE))

# replace NAs with blanks
rTrans[is.na(rTrans)] <- " "
rTrans <- as.data.frame(lapply(rTrans, as.numeric))

# Determine Response Model (RDA vs CCA)
decorana(rTrans)



#====================================================================
# Working Area #

