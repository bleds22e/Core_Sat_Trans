# Transient & NDVI experimentation
# EKB
# Sept 2017

# LIBRARIES #

library(dplyr)
library(ggplot2)
library(zoo)
library(gridExtra)
library(stringr)
source('scripts/data_processing.R')
source('scripts/RodentAbundances.R')

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

rdat <- read.csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Git/PortalData/Rodents/Portal_rodent.csv", na.strings = '')
NDVI <- read.csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Git/PortalData/NDVI/monthly_NDVI.csv")

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

par(mfrow = c(3,1))
hist(comm1_persist$rel_reg_persist)
hist(comm2_persist$rel_reg_persist)
hist(comm3_persist$rel_reg_persist)

#dev.copy(png, "plots/hist_portal_comms.png")
#dev.off()

# transients for each community

transients1 <- comm1_persist %>% 
  filter(rel_reg_persist < 1/3, species != "PB") %>% 
  select(species)
transients1 <- unlist(apply(transients1, 1, list), recursive = FALSE)

transients2 <- comm2_persist %>% 
  filter(rel_reg_persist < 1/3) %>% 
  select(species)
transients2 <- unlist(apply(transients2, 1, list), recursive = FALSE)

transients3 <- comm3_persist %>% 
  filter(rel_reg_persist < 1/3) %>% 
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


NDVI1 <- filter(NDVI, date %in% unique(abund_dates1$date))
NDVI2 <- filter(NDVI, date %in% unique(abund_dates2$date))
NDVI3 <- filter(NDVI, date %in% unique(abund_dates3$date))


# PLOT DATA #

NDVI1$date = as.yearmon(NDVI1$date)
NDVI2$date = as.yearmon(NDVI2$date)
NDVI3$date = as.yearmon(NDVI3$date)
abund_dates1$date = as.yearmon(abund_dates1$date)
abund_dates2$date = as.yearmon(abund_dates2$date)
abund_dates3$date = as.yearmon(abund_dates3$date)

(comm2_plot <- ggplot(data = NDVI2, aes(x = date, y = ndvi)) +
  geom_line(aes(color = 'NDVI')) +
  geom_point(data = abund_dates2, aes(x = date, y = transient_rel_abund, color = 'rel_abund'), size = 1) +
  ggtitle('2000-2009') +
  scale_color_manual(values = c(NDVI = '#333333', rel_abund = '#00CC00')) +
  ylab('NDVI value / Rel.Abund') +
  theme_bw())

#transient_comms <- grid.arrange(comm1_plot, comm2_plot, comm3_plot, nrow = 3)
#ggsave("plots/transient_comms_noPBcomm1.png", plot = transient_comms)

#===========================================================

# NDVI LENGTH # 

library(forecast)

NDVI_ts = ts(NDVI$ndvi, start = c(1992, 3), end = c(2015, 12), frequency = 12)
NDVI1_ts = ts(NDVI1$ndvi, start = c(1992, 3), end = c(1999, 12), frequency = 12)
NDVI2_ts = ts(NDVI2$ndvi, start = c(2000, 1), end = c(2009, 12), frequency = 12)
NDVI3_ts = ts(NDVI3$ndvi, start = c(2010, 1), end = c(2015, 12), frequency = 12)
comm1_ts = ts(abund_dates1$transient_rel_abund, start = c(1984, 2), end = c(1999, 12), frequency = 12)
comm2_ts = ts(abund_dates2$transient_rel_abund, start = c(2000, 1), end = c(2009, 12), frequency = 12)
comm3_ts = ts(abund_dates3$transient_rel_abund, start = c(2010, 1), end = c(2015, 12), frequency = 12)

NDVI_peak <- mutate(NDVI, NDVIpeak = ndvi - median(ndvi))
NDVI_peak$date = as.yearmon(NDVI_peak$date)

NDVI1_peak <- filter(NDVI_peak, date %in% unique(abund_dates1$date))
NDVI1_peak_ts = ts(NDVI1_peak$NDVIpeak, start = c(1992, 3), end = c(1999, 12), frequency = 12)
NDVI2_peak <- filter(NDVI_peak, date %in% unique(abund_dates2$date))
NDVI2_peak_ts = ts(NDVI2_peak$NDVIpeak, start = c(2000, 1), end = c(2009, 12), frequency = 12)
NDVI3_peak <- filter(NDVI_peak, date %in% unique(abund_dates3$date))
NDVI3_peak_ts = ts(NDVI3_peak$NDVIpeak, start = c(2010, 1), end = c(2015, 12), frequency = 12)

NDVI_all_peak <- bind_rows(NDVI1_peak, NDVI2_peak, NDVI3_peak)
NDVI_all_peak$date = as.yearmon(NDVI_all_peak$date)

par(mfrow = c(3,2))
ccf(NDVI1_peak_ts, comm1_ts)
plot(NDVI1_peak_ts)
lines(comm1_ts, col = 'blue')
points(comm1_ts, col = 'blue')
ccf(NDVI2_ts, comm2_ts)
plot(NDVI2_peak_ts)
lines(comm2_ts, col = 'blue')
points(comm2_ts, col = 'blue')
ccf(NDVI3_ts, comm3_ts)
plot(NDVI3_peak_ts)
par(new = TRUE)
lines(comm3_ts, col = 'blue')
points(comm3_ts, col = 'blue')

#dev.copy(png, "plots/ccf_and_relabund_plots.png")
#dev.off()

# plots

# NDVI1_peak$date = as.yearmon(NDVI1_peak$date)
# NDVI2_peak$date = as.yearmon(NDVI2_peak$date)
# NDVI3_peak$date = as.yearmon(NDVI3_peak$date)
# 
# (comm2_plot <- ggplot(data = NDVI2_peak, aes(x = date, y = NDVIpeak)) +
#     geom_line(aes(color = 'NDVIpeak')) +
#     geom_point(data = abund_dates2, aes(x = date, y = transient_rel_abund, color = 'rel_abund'), size = 1) +
#     ggtitle('2000-2009') +
#     scale_color_manual(values = c(NDVIpeak = '#333333', rel_abund = '#00CC00')) +
#     ylab('NDVI value / Rel.Abund') +
#     theme_bw())
# 
# transient_comms <- grid.arrange(comm1_plot, comm2_plot, comm3_plot, nrow = 3)
# ggsave("plots/transient_comms_noPBcomm1.png", plot = transient_comms)


