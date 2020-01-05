weather <- portalr::weather(level = "monthly", fill = TRUE)
weather_lunar <- portalr::weather(level = "newmoon", fill = TRUE)
ndvi <- portalr::ndvi(level = "monthly", fill = TRUE)
ndvi_lunar <- portalr::ndvi(level = "newmoon", fill = TRUE)

ggplot(ndvi, aes(x = date, y = ndvi)) +
  geom_line()

ggplot(ndvi_lunar) +
  geom_line(aes(x = newmoonnumber, y = ndvi)) 

ggplot(weather_lunar) +
  geom_line(aes(x = date, y = precipitation))

abund_lunar <- abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                                           plots = "all", unknowns = F, shape = "flat", 
                                           time = c("newmoon"), min_plots = 10)

moon_dates <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"))

# get rodent abundance data
abund <- portalr::abundance(path = "repo", level = "Site", type = "Rodents",
                            plots = "all", unknowns = F, shape = "flat", 
                            time = c("newmoon"), min_plots = 10)

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

trapping_lunar <- full_join(moon_dates, trapping_effort)

# add month, year, and effor to relative transient abundance
rel_abund_data <- full_join(trapping_lunar, rel_abund_trans)

# double the number of transients if approriate
for (i in 1:nrow(rel_abund_data)){
  if (rel_abund_data$nplots[i] <= 17){
    rel_abund_data$trans_rel_abund[i] <- rel_abund_data$trans_rel_abund[i]*2
  }
}

#====================================================
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
