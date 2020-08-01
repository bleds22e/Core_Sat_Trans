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

ndvi <- portalr::ndvi(level = "newmoon", fill = TRUE) %>% 
  filter(newmoonnumber < 452) # filter out forecasted ndvi data--real data ends at 2013
abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                            plots = "all", unknowns = F, shape = "flat", 
                            time = c("newmoon"), min_plots = 10)

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
  print(transients)
  
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

#===============================================================
# GET TRANSIENT RODENT ABUNDANCES 

# pre 1984, 1984-2000, 2000-2009, post-2009 - USING MOON DATES
comm1 <- filter(abund, newmoonnumber <= 80) # pre-1984
comm2 <- filter(abund, newmoonnumber > 80 & newmoonnumber < 279) # 1984-2000
comm3 <- filter(abund, newmoonnumber >= 279 & newmoonnumber < 403) # 2000-2009
comm4 <- filter(abund, newmoonnumber >= 403) # post-2009

# get transient relative abundance
rel_abund_trans1 <- get_transient_rel_abund(comm1)
rel_abund_trans2 <- get_transient_rel_abund(comm2)
rel_abund_trans3 <- get_transient_rel_abund(comm3)
rel_abund_trans4 <- get_transient_rel_abund(comm4)

rel_abund_trans <- do.call("rbind", list(rel_abund_trans1,
                                         rel_abund_trans2,
                                         rel_abund_trans3,
                                         rel_abund_trans4))

# get number of plots trapped per period
trapping_effort <- trapping %>% 
  select(period, sampled) %>% 
  group_by(period) %>% 
  summarise(nplots = sum(sampled)) 

# merge trapping effort, new moon dates, and transient rel. abundnace
trapping_lunar <- full_join(moon_dates, trapping_effort)
rel_abund_data <- full_join(trapping_lunar, rel_abund_trans)

# double the number of transients if appropriate
for (i in 1:nrow(rel_abund_data)){
  if (!is.na(rel_abund_data$nplots[i])){
    if (rel_abund_data$nplots[i] <= 17){
      rel_abund_data$trans_rel_abund[i] <- rel_abund_data$trans_rel_abund[i]*2
    } 
  }
}

#=========================================================
# COMBINE NDVI AND TRANSIENT DATA

# add transient relative abundance
NDVI_transient <- full_join(rel_abund_data, ndvi)

# fill in missing transient values
for (i in 2:nrow(NDVI_transient)){
  if (is.na(NDVI_transient$trans_rel_abund[i])){
    NDVI_transient$trans_rel_abund[i] <- (NDVI_transient$trans_rel_abund[i-1] + NDVI_transient$trans_rel_abund[i+1])/2
  }
}

# make function for this later
trans_value <- (NDVI_transient$trans_rel_abund[108] - NDVI_transient$trans_rel_abund[105])/3
NDVI_transient$trans_rel_abund[106] <- NDVI_transient$trans_rel_abund[105] + trans_value
NDVI_transient$trans_rel_abund[107] <- NDVI_transient$trans_rel_abund[105] + (2*trans_value)

NDVI_transient$trans_rel_abund[273:274] <- 0

#==========================================================
# MAKE NEW DATAFRAME FOR HEATMAP #

# make outline of data frame
nrows <- length(NDVI_transient$newmoonnumber)
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

for (i in 2:length(NDVI_transient$newmoonnumber)){
  
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
#ggsave("plots/newmoons_plots/heatmap_points.png")

ggplot(data = lags_long, aes(x = NDVI, y = d.NDVI, z = transients)) + 
  stat_summary_2d() +
  geom_vline(xintercept = median_ndvi, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.18)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/newmoons_plots/heatmap_summary.png")

ggplot(data = lags_long, aes(x = round(NDVI,2), y = round(d.NDVI, 2))) + 
  geom_tile(aes(fill = round(transients, 2))) + 
  geom_vline(xintercept = median_ndvi, col = "gray") +
  geom_hline(yintercept = 0, col = "gray") +
  scale_fill_viridis_c(limits = c(0, 0.18)) +
  facet_wrap(. ~ time_lag, nrow = 3, ncol = 4) +
  theme_bw()
#ggsave("plots/newmoons_plots/heatmap_2digit.png")

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
#ggsave("plots/newmoons_plots/transients_through_time_by_quadrants_0.05.png")

### Get Timeseries Plots ###

# pre 1984, 1984-2000, 2000-2009, post-2009 - USING MOON DATES
# comm1 <- filter(abund, newmoonnumber <= 80) # pre-1984
# comm2 <- filter(abund, newmoonnumber > 80 & newmoonnumber < 279) # 1984-2000
# comm3 <- filter(abund, newmoonnumber >= 279 & newmoonnumber < 403) # 2000-2009
# comm4 <- filter(abund, newmoonnumber >= 403) # post-2009

NDVI_transient$community <- NA
for (i in 1:length(NDVI_transient$community)) {
  
  if (is.na(NDVI_transient$newmoonnumber[i])) {
    NDVI_transient$community[i] = NA
  } else if (NDVI_transient$newmoonnumber[i] <= 80) {
    NDVI_transient$community[i] = 'comm1'
  } else if (NDVI_transient$newmoonnumber[i] > 80 & NDVI_transient$newmoonnumber[i] < 279) {
    NDVI_transient$community[i] = "comm2"
  } else if (NDVI_transient$newmoonnumber[i] >= 279 & NDVI_transient$newmoonnumber[i] < 403) {
    NDVI_transient$community[i] = "comm3"
  } else if (NDVI_transient$newmoonnumber[i] >= 403) {
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
NDVI_transient$ndvi_median <- NA
for (i in 1:length(NDVI_transient$ndvi_median)) {
  if (!is.na(NDVI_transient$ndvi[i])){
    NDVI_transient$ndvi_median[i] <- NDVI_transient$ndvi[i] - median
  } else {
    NDVI_transient$ndvi_median[i] <- NA
  }
}

NDVI_transient$newmoondate <- as.Date(NDVI_transient$newmoondate)

ggplot(NDVI_transient) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_line(aes(x = newmoondate, y = ndvi_median, color = "NDVI")) +
  geom_line(aes(x = newmoondate, y = trans_rel_abund, color = "Transients")) +
  facet_wrap(vars(community), scales = "free", nrow = 4) +
  scale_color_manual("",
                     breaks = c("NDVI", "Transients"),
                     values = c("gray", "black")) +
  ylab('NDVI value / Rel.Abund') +
  xlab('Date') +
  theme_bw()
#ggsave("plots/newmoons_plots/timeseries_by_community.png")

#================================================================
# RUN EDM and CCF -- script from transient_convergent_crossmapping.R

new_df <- NDVI_transient[-c(1:49, 452:523),] %>% 
  ungroup() %>% 
  select(newmoondate, trans_rel_abund, ndvi)

# CCF

NDVI <- new_df$ndvi
Transient_Rel_Abund <- new_df$trans_rel_abund

acf(NDVI)
lag.plot(new_df$ndvi)
pacf(new_df$ndvi)
fit <- forecast::auto.arima(new_df$ndvi)

acf(Transient_Rel_Abund)
lag.plot(new_df$trans_rel_abund)
pacf(new_df$trans_rel_abund)
fit2 <- forecast::auto.arima(new_df$trans_rel_abund)

par(mfrow = c(3,1))
acf(NDVI)
acf(Transient_Rel_Abund)
ccf(NDVI, Transient_Rel_Abund)

dev.copy(png, "plots/acf_ccf.png")
dev.off()

# CONVERGENT CROSS-MAPPING
library(rEDM)

ts <- new_df$trans_rel_abund
lib <- c(1, length(ts))
pred <- c(1, length(ts))
simplex_output <- simplex(ts, lib, pred, E = 1:15, silent = TRUE)
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")


smap_output <- list(s_map(ts, lib, pred, E = 11, silent = TRUE))
plot(smap_output[[1]]$theta, smap_output[[1]]$rho, type = "l", xlim = c(0, 4), 
     xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

vars <- colnames(new_df[2:3])
n <- nrow(new_df)

ccm_rho_matrix <- matrix(NA, nrow = length(vars), ncol = length(vars), dimnames = list(vars,vars))

for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    out_temp <- ccm(new_df, E = 6, lib_column = ccm_from, target_column = ccm_to, 
                    lib_sizes = n, replace = FALSE, silent = TRUE)
    ccm_rho_matrix[ccm_from, ccm_to] <- out_temp$rho
  }
}


corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars))

for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    cf_temp <- ccf(new_df[, ccm_from], new_df[, ccm_to], type = "correlation", 
                   lag.max = 12, plot = FALSE)$acf
    corr_matrix[ccm_from, ccm_to] <- max(abs(cf_temp))
  }
}

head(ccm_rho_matrix)
head(corr_matrix)

trans_xmap_ndvi <- ccm(new_df, E = 6, random_libs = TRUE, lib_column = "trans_rel_abund", 
                       target_column = "ndvi", lib_sizes = seq(10, 400, by = 10), num_samples = 500, 
                       silent = TRUE)
ndvi_xmap_trans <- ccm(new_df, E = 6, random_libs = TRUE, lib_column = "ndvi", 
                       target_column = "trans_rel_abund", lib_sizes = seq(10, 400, by = 10), num_samples = 500, 
                       silent = TRUE)
ccm_out <- list(ccm_means(trans_xmap_ndvi), ccm_means(ndvi_xmap_trans))

plot(ccm_out[[1]]$lib_size, pmax(0, ccm_out[[1]]$rho), type = "l", col = "red",  
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(ccm_out[[2]]$lib_size, pmax(0, ccm_out[[2]]$rho), col = "blue")
abline(h = corr_matrix['trans_rel_abund', 'ndvi'], col = "black", lty = 2)
legend(x = "topright", legend = c("trans_xmap_ndvi", "ndvi_xmap_trans"),
       col = c("red", "blue"), lwd = 1, inset = 0.02)

params <- expand.grid(lib_column = vars, target_column = vars, tp = -10:15)
params <- params[params$lib_column != params$target_column, ]
params$E <- 6

output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
  ccm(new_df, E = params$E[i], lib_sizes = NROW(new_df), 
      random_libs = FALSE, lib_column = params$lib_column[i], target_column = params$target_column[i], 
      tp = params$tp[i], silent = TRUE)
}))

output$direction <- paste(output$lib_column, "xmap to\n", output$target_column)

time_delay_ccm_fig <- ggplot(output, aes(x = tp, y = rho, color = direction)) + 
  geom_line() + theme_bw()
print(time_delay_ccm_fig)

