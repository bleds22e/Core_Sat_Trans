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
trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))

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
  select(species)
transients1 <- unlist(apply(transients1, 1, list), recursive = FALSE)

transients2 <- comm2_persist %>% 
  filter(rel_reg_persist <= 1/3 ) %>% 
  select(species)
transients2 <- unlist(apply(transients2, 1, list), recursive = FALSE)

transients3 <- comm3_persist %>% 
  filter(rel_reg_persist <= 1/3) %>% 
  select(species)
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

dates <- select(rdat, year, month, period) %>% 
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
  select(period, plot, sampled) %>% 
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
NDVI <- arrange(NDVI, year, month)

for (i in 1:nrow(NDVI)){
  if (str_length(NDVI$month[i]) == 1){
    NDVI$month[i] = paste('0', NDVI$month[i], sep = '')
  }
}

NDVI <- select(NDVI, year, month, ndvi) %>% 
  tidyr::unite(date, year, month, sep = '-') %>% 
  distinct()

# make a column for NDVI centered around the median
NDVI_peak <- mutate(NDVI, NDVIpeak = ndvi - median(ndvi))

# filter NDVI by community dates
NDVI1 <- filter(NDVI_peak, date %in% unique(abund_dates1$date))
NDVI2 <- filter(NDVI_peak, date %in% unique(abund_dates2$date))
NDVI3 <- filter(NDVI_peak, date %in% unique(abund_dates3$date))
NDVI_all <- bind_rows(NDVI1, NDVI2, NDVI3)


### COMBINE TRANSIENTS AND NDVI ###

combined_data <- full_join(abund_dates, NDVI_all, by = "date")
combined_data$date = as.yearmon(combined_data$date)

# deal with periods that have two different months
# select first month out of the two?

combined_data <- combined_data[!duplicated(combined_data$period),]

# get all yearmon dates
all_dates <- as.data.frame(seq(as.yearmon(min(NDVI_peak$date)), as.yearmon(max(NDVI_peak$date)), by = 1/12))
colnames(all_dates) = c("date")
additional_year <- as.data.frame(seq(from = max(all_dates$date) + 1/12, by = 1/12, length.out = 12))
colnames(additional_year) = c("date")
all_dates <- rbind(all_dates, additional_year)

combined_data <- full_join(combined_data, all_dates, by = "date") %>% 
  arrange(date)
combined_data <- combined_data[-c(1:65),]

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


# INTERPOLATE MISSING VALUES #

for (i in 1:nrow(combined_data)){

  # fill in missing transient values
  if (is.na(combined_data$transient_rel_abund[i])){
    combined_data$transient_rel_abund[i] <- (combined_data$transient_rel_abund[i-1] + combined_data$transient_rel_abund[i+1])/2
  }
  
  # fill in missing NDVI values
  if (is.na(combined_data$NDVIpeak[i])){
    combined_data$NDVIpeak[i] <- (combined_data$NDVIpeak[i-1] + combined_data$NDVIpeak[i+1])/2
  } 
  
}

combined_data$transient_rel_abund[165:166] <- 0
combined_data$NDVIpeak[165:166] <- (combined_data$NDVIpeak[164] + combined_data$NDVIpeak[167])/2


# FIND NDVI PEAKS #

# quantiles
upper <- as.numeric(quantile(NDVI_peak$NDVIpeak, .75))
lower <- as.numeric(quantile(NDVI_peak$NDVIpeak, .25))

combined_data$NDVI_quantiles = NA

for (i in 1:nrow(combined_data)){
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

for (i in 1:nrow(combined_data)){
  
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
                     pulseDuration = numeric(nrows),
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

pulse_seq <- as.numeric(seq(from = 1, to = nrows))
pulses$pulseID <- pulse_seq

for (i in 1:length(pulse_seq)){
  
  # select only rows from that pulse
  current_pulse <- combined_pulses[combined_pulses$pulseID == i,]
  
  # find the maximum NDVI value
  pulses$pulseMax[i] <- max(current_pulse$NDVIpeak)
  
  # find pulse duration
  pulses$pulseDuration[i] <- nrow(current_pulse)
  
  # find pulse sum
  pulses$pulseSum[i] <- sum(current_pulse$NDVIpeak)
  
  # find time (months) between end of current pulse and start of next pulse
  pulses$timeBetween[i] <- 12 * (min(combined_pulses$date[combined_pulses$pulseID == i + 1]) - max(current_pulse$date))
  
  # find the transient rel abundances
  min_date <- min(current_pulse$date)
  pulses$transientT0[i] <- combined_data$transient_rel_abund[combined_data$date == min_date]
  pulses$transientT1[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 1/12)]
  pulses$transientT2[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 2/12)]
  pulses$transientT3[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 3/12)]
  pulses$transientT4[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 4/12)]
  pulses$transientT5[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 5/12)]
  pulses$transientT6[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 6/12)]
  pulses$transientT7[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 7/12)]
  pulses$transientT8[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 8/12)]
  pulses$transientT9[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 9/12)]
  pulses$transientT10[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 10/12)]
  pulses$transientT11[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 11/12)]
  pulses$transientT12[i] <- combined_data$transient_rel_abund[combined_data$date == (min_date + 12/12)]
  
}


#========================================================================
# Run CCA Analysis

# separate into env and response variables
pulses_NDVI <- pulses[,1:5]
pulses_transients <- pulses[,6:18]

# remove rows where rowSum is 0
rsum <- rowSums(pulses_transients)
pulses_transients <- pulses_transients[-c(29,30,33),]

# histograms of transient distributions
mapply(hist, as.data.frame(pulses_transients[,1:13], 
                           main = colnames(pulses_transients[,1:13]),
                           xlab = "abundance"))

# transformations
log.full <- log1p(pulses_transients)
hell.full <- decostand(pulses_transients, "hellinger")


# check row and column sum variability
rsum <- rowSums(hell.full, na.rm = TRUE)
csum <- colSums(hell.full, na.rm = TRUE)
hist(rsum)
hist(csum)
cv(rsum)
cv(csum)

# Determine Response Model (RDA vs CCA)
decorana(rTrans)
  # need to run RDA instead of CCA (after change to 1/3 rather than 1/2)

# scale the explanatory variables
vars <- pulses_NDVI[-c(29,30,33), -1]
vars$timeBetween <- log(vars$timeBetween) 
cv(colSums(vars))
vars_z <- scores(vars)

explanatory <- as.data.frame(scale(vars_z))
mapply(hist, as.data.frame(vars, 
                           main = colnames(vars)))


# NOTES FOR RUNNING THE CCA
#   do the pulse variables need to be transformed? -- maybe log transform 1 and 4?
#   probably need to z-standardize them anyway (check other labs' code)
#   check for correlation in the variables (pairwise)

# trying RDA
trans.rda <- rda(hell.full ~ ., explanatory)
summary(trans.rda)

R2 <- RsquareAdj(trans.rda)$r.squared
R2adj <- RsquareAdj(trans.rda)$adj.r.squared

############### NEED TO SWITCH TO RDA ######################
# run CCA

ca <- cca(rTrans)
plot(ca)
summary(ca)

transient.CCA <- cca(rTrans ~ ., data = explanatory)
summary(transient.CCA)

anova.cca(transient.CCA)
anova.cca(transient.CCA, by = 'axis')

par(mfrow = c(1,2))
plot(transient.CCA$CCA$wa[,1], transient.CCA$CCA$wa[,2], 
     xlab = "CCA Axis 1", ylab = "CCA Axis 2")
plot(transient.CCA$CCA$u[,1], transient.CCA$CCA$u[,2], 
     xlab = "CCA Axis 1", ylab = "CCA Axis 2")

spenvcor(transient.CCA)

transient.CCA$CCA$biplot

plot(transient.CCA, choices=c(1,2), display = c('wa','sp','bp'), scaling = 2)
plot(transient.CCA, choices=c(1,2), display = c('lc','sp','bp'), scaling = 2)

#====================================================================
# Working Area #

