# Spordic species abundance ~ Monthly NDVI values
# May 1, 2017
# EKB

# maybe try by year or season? 


# LIBRARIES and SOURCE CODE #

library(stringr)
library(ggplot2)
library(zoo)
library(gridExtra)
source('scripts/data_processing.R')
source('scripts/RodentAbundances.R')

# GET DATA #

rdat <- read.csv("~/Dropbox (UFL)/Git/PortalData/Rodents/Portal_rodent.csv", na.strings = '')
NDVI <- read.csv("~/Dropbox (UFL)/Git/PortalData/NDVI/monthly_NDVI.csv")

#=======================================================
# BY PERIOD #

# All Sporadic Summed Together #

abund <- abundance(path = "repo", level = "Site", type = "Rodents",
                   length = "all", unknowns = F, incomplete = F,
                   shape = "flat", time = "period")

total_by_period <- abund %>% group_by(period) %>% 
                   summarise(total = sum(abundance))

# different sporadic groupings
transients <- c('BA', 'PH', 'PI', 'PL', 'PM', 'RF', 'RO', 'SF', 'SH', 'SO')
transients_PF <- c('BA', 'PF', 'PH', 'PI', 'PL', 'PM', 'RF', 'RO', 'SF', 'SH', 'SO')
transients_RM <- c('BA', 'PH', 'PI', 'PL', 'PM', 'RM', 'RF', 'RO', 'SF', 'SH', 'SO')
transients_PF_RM <- c('BA', 'PF', 'PH', 'PI', 'PL', 'PM', 'RM', 'RF', 'RO', 'SF', 'SH', 'SO')

# get relative abundances
abund <- filter(abund, species %in% transients) %>% 
         group_by(period) %>% 
         summarise(transients = sum(abundance))

abund <- inner_join(abund, total_by_period)
rel_abund_transients <- mutate(abund, transient_rel_abund = transients/total)

# add date to abundances
for (i in 1:length(rdat$mo)){
  if (str_length(rdat$mo[i]) == 1){
    rdat$mo[i] = paste('0', rdat$mo[i], sep = '')
  }
}

dates <- select(rdat, yr, mo, period) %>% 
  unite(Date, yr, mo, sep = '-') %>% 
  distinct()

# filter abundances by available NDVI
NDVI = NDVI[-(1:110),-1]
abund_dates <- inner_join(rel_abund_transients, dates, by = 'period') %>% 
  filter(Date %in% unique(NDVI$Date))

# plot data 

NDVI$Date = as.yearmon(NDVI$Date)
abund_dates$Date = as.yearmon(abund_dates$Date)

(spor <- ggplot(data = NDVI, aes(x = Date, y = NDVI)) +
  geom_line(aes(color = 'NDVI')) +
  geom_point(data = abund_dates, aes(x = Date, y = transient_rel_abund, color = 'Rel.Abund'), size = 1) +
  scale_color_manual(values = c(NDVI = '#333333', Rel.Abund = '#00CC00'))+
  ggtitle('Sporadics') +
  ylab('NDVI value / Rel.Abund')+
  theme_bw())

four_spor_levels <- grid.arrange(spor, spor_PF, spor_RM, spor_PF_RM, nrow = 4)
#ggsave("plots/four_sporadic_levels_NDVI.png", plot = four_spor_levels)

# sporadics ~ NDVI

NDVI_abund <- full_join(abund_dates, NDVI, by = 'Date')
ggplot(data = NDVI_abund, aes(x = NDVI, y = transient_rel_abund)) +
  geom_point() +
  theme_bw()
#ggsave("plots/NDVI_vs_transientRelAbund.png")

##########

# Each Transient Species Separately #
abund <- abundance(path = "repo", level = "Site", type = "Rodents",
                   length = "all", unknowns = F, incomplete = F,
                   shape = "flat", time = "period")

total_by_period_each <- abund %>% group_by(period) %>% 
  summarise(total = sum(abundance))

transients <- c('BA', 'PF', 'PH', 'PI', 'PL', 'PM', 'RF', 'RM', 'RO', 'SF', 'SH', 'SO')
abund <- filter(abund, species %in% transients) %>% 
  group_by(period)

abund <- left_join(abund, total_by_period)
rel_abund_transients <- mutate(abund, rel_abund = abundance/total)

abund_dates <- left_join(rel_abund_transients, dates, by = 'period')
abund_dates$Date = as.yearmon(abund_dates$Date)
abund_dates <- filter(abund_dates, Date %in% unique(NDVI$Date))

NDVI_abund <- full_join(abund_dates, NDVI, by = 'Date') %>% 
              filter(species != 'NA')

ggplot(data = NDVI_abund, aes(x = Date, y = NDVI)) +
    geom_line(aes(color = 'NDVI')) +
    geom_point(aes(x = Date, y = rel_abund, color = 'Rel.Abund'), size = 1) +
    facet_wrap(~species, nrow = 4, ncol = 3)+
    scale_color_manual(values = c(NDVI = '#333333', Rel.Abund = '#00CC00'))+
    ylab('NDVI value / Rel.Abund')+
    theme_bw()
#ggsave("plots/NDVI_sporadics_by_sp.png")

#=======================================================
# BY SEASON #

abund <- abundance(path = "repo", level = "Site", type = "Rodents",
                   length = "all", unknowns = F, incomplete = F,
                   shape = "flat", time = "period")

season_dates <- select(rdat, month = mo, year = yr, period) %>% 
                distinct()
abund_season_dates <- left_join(abund, season_dates)

abund_seasonal <- add_seasons(abund_season_dates, 2)

# get seasonal relative abundances

total_by_season <- abund_seasonal %>% group_by(season2_yr, season2) %>% 
  summarise(total = sum(abundance))

transients <- c('BA', 'PF', 'PH', 'PI', 'PL', 'PM', 'RF', 'RM', 'RO', 'SF', 'SH', 'SO')

rel_abundances <- filter(abund_seasonal, species %in% transients) %>% 
  group_by(season2_yr, season2) %>% 
  summarise(transients = sum(abundance))

seasonal_abund <- inner_join(rel_abundances, total_by_season)
seasonal_rel_abund <- mutate(seasonal_abund, transient_rel_abund = transients/total)

# average NDVI by season

NDVI <- read.csv("C:/Users/ellen.bledsoe/Desktop/Git/PortalData/NDVI/monthly_NDVI.csv")
NDVI = NDVI[-(1:110),-1]
NDVI <- tidyr::separate(NDVI, col = Date, into = c("year", "month"), sep = '-')
NDVI$year <- as.numeric(NDVI$year)
NDVI$month <- as.numeric(NDVI$month)

NDVI_season <- add_seasons(NDVI, 2)
