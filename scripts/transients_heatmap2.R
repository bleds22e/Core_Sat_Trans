# Transients Heatmaps 2
# EKB
# February 2019

# LIBRARIES and DATA #

library(RCurl)
library(zoo)
library(tidyverse)
library(portalr)

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

#=========================================================
# COMBINE NDVI AND TRANSIENT DATA

# add transient relative abundance
NDVI_transient <- full_join(rel_abund_data, NDVI)

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

# make function for this later
trans_value <- (NDVI_transient$trans_rel_abund[107] - NDVI_transient$trans_rel_abund[104])/3
NDVI_transient$trans_rel_abund[105] <- NDVI_transient$trans_rel_abund[104] + trans_value
NDVI_transient$trans_rel_abund[106] <- NDVI_transient$trans_rel_abund[104] + (2*trans_value)

NDVI_transient$trans_rel_abund[267:268] <- 0

#==========================================================
# MAKE NEW DATAFRAME FOR HEATMAP #

NDVI_transient$date <- as.yearmon(paste(NDVI_transient$year, NDVI_transient$month), "%Y %m") 

# make outline of data frame
nrows <- length(NDVI_transient$date)
lags <- data.frame(NDVI = NDVI_transient$ndvi,
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

for (i in 2:length(NDVI_transient$date)){
  
  # find change in NDVI
  lags$d.NDVI[i] <- lags$NDVI[i] - lags$NDVI[i-1]
  
  # find the transient rel abundances
  lags$transientT0[i] <- NDVI_transient$trans_rel_abund[i]
  lags$transientT1[i] <- NDVI_transient$trans_rel_abund[i + 1]
  lags$transientT2[i] <- NDVI_transient$trans_rel_abund[i + 2]
  lags$transientT3[i] <- NDVI_transient$trans_rel_abund[i + 3]
  lags$transientT4[i] <- NDVI_transient$trans_rel_abund[i + 4]
  lags$transientT5[i] <- NDVI_transient$trans_rel_abund[i + 5]
  lags$transientT6[i] <- NDVI_transient$trans_rel_abund[i + 6]
  lags$transientT7[i] <- NDVI_transient$trans_rel_abund[i + 7]
  lags$transientT8[i] <- NDVI_transient$trans_rel_abund[i + 8]
  lags$transientT9[i] <- NDVI_transient$trans_rel_abund[i + 9]
  lags$transientT10[i] <- NDVI_transient$trans_rel_abund[i + 10]
  lags$transientT11[i] <- NDVI_transient$trans_rel_abund[i + 11]
  
}

lags <- lags[complete.cases(lags),]

#=================================================================
# PLOTTING #

median_ndvi <- median(NDVI_transient$ndvi, na.rm = TRUE)

### Heatmap style plotting ###

# rearrange data for easier plotting
lags_long <- tidyr::gather(lags, "time_lag", "transients", 3:14)

new_levels = c("transientT0", "transientT1", "transientT2",
               "transientT3", "transientT4", "transientT5",
               "transientT6", "transientT7", "transientT8",
               "transientT9", "transientT10", "transientT11")
lags_long <- arrange(mutate(lags_long, time_lag = factor(time_lag, levels = new_levels)), time_lag)

# various heatmap plots
ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  stat_summary_2d() +
  geom_point(shape = 1, col = "white") + 
  geom_vline(xintercept = median_ndvi, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.18)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/GIMMs_plots/heatmap_points.png")

ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  stat_summary_2d() +
  geom_vline(xintercept = median_ndvi, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.18)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/GIMMs_plots/heatmap_summary.png")

ggplot(data = lags_long, aes(x = round(NDVI,2), y = round(d.NDVI, 2))) + 
  geom_tile(aes(fill = round(transients, 2))) + 
  geom_vline(xintercept = median_ndvi, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.18)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/GIMMs_plots/heatmap_2digit.png")

### Get "contingency" tables ###

# get quadrants
x.neg_y.pos <- lags_long %>% 
  filter(NDVI < median_ndvi, d.NDVI > 0, transients > 0.05) %>% 
  group_by(time_lag) %>% 
  summarise(x.neg_y.pos = mean(transients))
x.pos_y.pos <- lags_long %>% 
  filter(NDVI >= median_ndvi, d.NDVI > 0, transients > 0.05) %>% 
  group_by(time_lag) %>% 
  summarise(x.pos_y.pos = mean(transients))
x.pos_y.neg <- lags_long %>% 
  filter(NDVI >= median_ndvi, d.NDVI <= 0, transients > 0.05) %>% 
  group_by(time_lag) %>% 
  summarise(x.pos_y.neg = mean(transients))
x.neg_y.neg <- lags_long %>% 
  filter(NDVI < median_ndvi, d.NDVI <= 0, transients > 0.05) %>% 
  group_by(time_lag) %>% 
  summarise(x.neg_y.neg = mean(transients))

quadrant_means <- plyr::join_all(list(x.neg_y.pos, x.pos_y.pos, x.pos_y.neg, x.neg_y.neg), by = "time_lag")
quadrant_means_long <- gather(quadrant_means, key = "quadrant", value = "mean_transients", 2:5)

ggplot(quadrant_means_long, aes(time_lag, mean_transients, color = quadrant, group = quadrant)) + 
  geom_point(size = 2) +
  geom_smooth() +
  theme_bw() +
  xlab("Lag Time") +
  ylab("Mean Transients (above 0.025)") +
  theme(axis.text.x = element_text(angle = -45, hjust = -.1))
#ggsave("plots/GIMMs_plots/transients_through_time_by_quadrants_0.05.png")

# library(ggpubr)
# plot_arranged <- ggarrange(plot_low, plot_mid, plot_high, nrow = 1, ncol = 3,
#           labels = c("transients > 0.01", "transients > 0.025", "transients > 0.05"))
# plot_arranged
#ggsave("plots/GIMMs_plots/transients_through_time_different_cutoffs.png")

### 3D plots ###
# library(spatstat)
# 
# point_pattern3 <- pp3(x = lags$NDVI, y = lags$d.NDVI, z = lags$transientT0, 
#                       as.box3(xrange = c(min(lags$NDVI), max(lags$NDVI)), 
#                               yrange = c(min(lags$d.NDVI), max(lags$d.NDVI)), 
#                               zrange = c(min(lags$transientT0), max(lags$transientT0))))
# plot.pp3(point_pattern3)
# plot.pp3(point_pattern3, legend = TRUE)
# 
# library(plotly)
# 
# plot <- plot_ly(lags, x = ~NDVI, y = ~d.NDVI, z = ~transientT6) %>% 
#   add_markers(fill = "tozeroy") %>% 
#   layout(scene = list(xaxis = list(title = "NDVI"),
#                       yaxis = list(title = "d.NDVI"),
#                       zaxis = list(title = "transients")))
# 
# plot_ly(x = lags$NDVI, y = lags$d.NDVI, z = lags$transientT5,
#         type = 'scatter3d', mode = "markers", color = lags$transientT5)

### Get Timeseries Plots ###

# add community designation to dataframe
comm1 <- filter(abund, period <= 75) # pre-1984
comm2 <- filter(abund, period > 75 & period < 261) # 1984-2000
comm3 <- filter(abund, period >= 261 & period < 381) # 2000-2009
comm4 <- filter(abund, period >= 381) # post-2009

NDVI_transient$community <- NA
for (i in 1:length(NDVI_transient$community)) {
  
  if (is.na(NDVI_transient$period[i])) {
    NDVI_transient$community[i] = NA
  } else if (NDVI_transient$period[i] <= 75) {
    NDVI_transient$community[i] = 'comm1'
  } else if (NDVI_transient$period[i] > 75 & NDVI_transient$period[i] < 261) {
    NDVI_transient$community[i] = "comm2"
  } else if (NDVI_transient$period[i] >= 261 & NDVI_transient$period[i] < 381) {
    NDVI_transient$community[i] = "comm3"
  } else if (NDVI_transient$period[i] >= 381) {
    NDVI_transient$community[i] = "comm4"
  }
  
  if (is.na(NDVI_transient$community[i])) {
    NDVI_transient$community[i] <- NDVI_transient$community[i-1]
  }
  
  if (is.na(NDVI_transient$community[i])) {
    NDVI_transient$community[i] <- NDVI_transient$community[i-1]
  }
  
}


median <- median(NDVI_transient$ndvi, na.rm = TRUE)
for (i in 1:length(NDVI_transient$ndvi_median)) {
  if (!is.na(NDVI_transient$ndvi[i])){
    NDVI_transient$ndvi_median[i] <- NDVI_transient$ndvi[i] - median
  } else {
    NDVI_transient$ndvi_median[i] <- NA
  }
}

ggplot(NDVI_transient) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_line(aes(x = date, y = ndvi_median, color = "NDVI")) +
  geom_line(aes(x = date, y = trans_rel_abund, color = "Transients")) +
  facet_wrap(vars(community), scales = "free", nrow = 4) +
  scale_color_manual("",
                     breaks = c("NDVI", "Transients"),
                     values = c("black", "green")) +
  ylab('NDVI value / Rel.Abund') +
  xlab('Date') +
  theme_bw()
#ggsave("plots/GIMMs_plots/timeseries_by_community.png")
