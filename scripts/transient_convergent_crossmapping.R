# Attempting Convergent Crossmapping
# EKB
# Oct. 2019
# first 180 lines copied from transients_heatmap2

library(zoo)
library(tidyverse)
library(RCurl)
library(rEDM)

rdat <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"), na.strings = '')
trapping <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"))
NDVI <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/NDVI/monthly_NDVI.csv"))

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

# add transient relative abundance -- landsat
# NDVI_transient <- full_join(rel_abund_data, NDVI) %>% 
#   rename("ndvi" = "x")

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

#=============================================================
# Attemping EDM and Convergent Crossmapping
#=============================================================

# landsat df
# new_df <- NDVI_transient[-c(1:178, 493:513),] %>% 
#   tibble::rowid_to_column() %>% 
#   ungroup() %>% 
#   select(rowid, trans_rel_abund, ndvi)

#gimms df
new_df <- NDVI_transient[-c(1:49, 452:523),] %>% 
  tibble::rowid_to_column() %>% 
  ungroup() %>% 
  select(rowid, trans_rel_abund, ndvi)

ts <- new_df$trans_rel_abund
lib <- c(1, length(ts))
pred <- c(1, length(ts))
simplex_output <- simplex(ts, lib, pred, E = 1:15, silent = TRUE)
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")


smap_output <- list(s_map(ts, lib, pred, E = 1, silent = TRUE))
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

# stopped here -- how do you determine library size?!
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
