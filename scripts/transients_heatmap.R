# transients heatmaps
# January 2019
# EKB

# Script is the same as transients_by_season.R until
# about row 216

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
#combined_data <- combined_data[-c(1:74),]
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


# MAKE NEW DATAFRAME FOR HEATMAP #

# make outline of data frame
nrows <- length(combined_data$date)
lags <- data.frame(NDVI = combined_data$NDVIpeak,
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

# plot heat maps

######
# try transforming by ways suggested in stack overflow
# and by centering around column (or maybe row) mean?
# also `rcompanion` transformTukey()
# bestNormalize

bN_NDVI <- bestNormalize(lags$NDVI)
bN_NDVI
oN_NDVI <- orderNorm(lags$NDVI)
par(mfrow = c(1,2))
hist(lags$NDVI)
hist(oN_NDVI$x.t)

bN_d.NDVI <- bestNormalize(lags$d.NDVI)
bN_d.NDVI
par(mfrow = c(1,2))
hist(lags$d.NDVI)
hist(bN_d.NDVI$x.t)

bN_transient <- bestNormalize(lags$transientT0)
bN_transient
oN_transient <- orderNorm(lags$transientT0)
par(mfrow = c(1,3))
hist(lags$transientT0)
hist(bN_transient$x.t)
hist(oN_transient$x.t)


# lags$NDVI_bin <- findInterval(lags$NDVI, seq(-0.1, 0.35, by = 0.025))
# lags$d.NDVI_bin <- findInterval(lags$d.NDVI, seq(-0.15, 0.375, by = 0.025))

lags_oN <- lags
for (i in 1:length(lags)) {
  i = i
  orderNorm_obj <- orderNorm(lags_oN[,i])
  lags_oN[,i] <- orderNorm_obj$x.t
}

par(mfrow = c(4,4))
for (i in 1:length(lags_oN)){
  i = i
  hist(lags_oN[,i])
}

# with no tranformation
lags_long <- tidyr::gather(lags, "time_lag", "transients", 3:14)

new_levels = c("transientT0", "transientT1", "transientT2",
               "transientT3", "transientT4", "transientT5",
               "transientT6", "transientT7", "transientT8",
               "transientT9", "transientT10", "transientT11")
lags_long <- arrange(mutate(lags_long, time_lag = factor(time_lag, levels = new_levels)), time_lag)

ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  stat_summary_2d() +
  geom_point(shape = 1, col = "white") + 
  geom_vline(xintercept = 0, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.35)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
ggsave("plots/heatmap_points.png")

ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  stat_summary_2d() +
  geom_vline(xintercept = 0, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.35)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
ggsave("plots/heatmap_summary.png")

ggplot(data = lags_long, aes(x = NDVI_bin, y = d.NDVI_bin)) + 
  geom_tile(aes(fill = transients)) + 
  geom_vline(xintercept = 4.5, col = "gray") +
  geom_hline(yintercept = 6.5, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.35)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
ggsave("plots/heatmap_geom_tile.png")

ggplot(data = lags_long, aes(x = round(NDVI,1), y = round(d.NDVI, 1))) + 
  geom_tile(aes(fill = round(transients, 1))) + 
  geom_vline(xintercept = 0, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.35)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
ggsave("plots/heatmap_1digit.png")


### with orderNorm tranformations

lags_oN$NDVI_bin <- findInterval(lags_oN$NDVI, seq(-3, 3, by = 0.5))
lags_oN$d.NDVI_bin <- findInterval(lags_oN$d.NDVI, seq(-3, 3, by = 0.5))


lags_long <- tidyr::gather(lags_oN, "time_lag", "transients", 3:14)

new_levels = c("transientT0", "transientT1", "transientT2",
               "transientT3", "transientT4", "transientT5",
               "transientT6", "transientT7", "transientT8",
               "transientT9", "transientT10", "transientT11")
lags_long <- arrange(mutate(lags_long, time_lag = factor(time_lag, levels = new_levels)), time_lag)

ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  stat_summary_2d() +
  geom_point(shape = 1, col = "white") + 
  geom_vline(xintercept = 0, col = "gray") +
  geom_hline(yintercept = .26, col = "gray") +
  scale_fill_viridis_c() +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/heatmap_points_oN.png")

ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  geom_vline(xintercept = 0, col = "gray") +
  geom_hline(yintercept = 0.26, col = "gray") +
  stat_summary_2d() +
  scale_fill_viridis_c() +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
ggsave("plots/heatmap_summary_oN.png")

ggplot(data = lags_long, aes(x = NDVI_bin, y = d.NDVI_bin)) + 
  geom_tile(aes(fill = transients)) + 
  geom_vline(xintercept = 6.5, col = "gray") +
  geom_hline(yintercept = 6.5, col = "gray") +
  scale_fill_viridis_c() +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
ggsave("plots/heatmap_geom_tile_oN.png")

ggplot(data = lags_long, aes(x = round(NDVI,0), y = round(d.NDVI, 0))) + 
  geom_tile(aes(fill = round(transients, 1))) + 
  geom_vline(xintercept = 0, col = "gray") +
  geom_hline(yintercept = 0.26, col = "gray") +
  scale_fill_viridis_c() +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/heatmap_1digit_oN.png")