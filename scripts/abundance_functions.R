#####################
# FUNCTIONS

# DATA functions

select_data <- function(data){
  # function for selecting only relevant columns
  dat <- select(data, year, plot, subplot, species)
  return(dat)
}

group_data <- function(data){
  # group and arrange a dataset by species, year, plot, and subplot
  dat <- data %>% 
    group_by(species, year, plot, subplot) %>% 
    arrange(species, year, plot, subplot)
  return(dat)
}

add_siteID <- function(data){
  # make a column with a unique location code
  dat <- data %>% mutate(siteID = paste(plot, subplot, sep = "_"))
  return(dat)
}

prep_data <- function(data){
  # run select_data, group_data, and add_siteID
  data1 <- select_data(data)
  data2 <- group_data(data1)
  data3 <- add_siteID(data2)
  return(data3)
}

adj_abundance <- function(data){
  # calculate the abundance per species, adjusted by the largest average
  abund <- data %>% select(year, siteID, species) %>% 
    group_by(species) %>% 
    summarise(count = n()) %>% 
    mutate(adj_abund = count/max(count))
  return(abund)
}

adj_avg_abund <- function(data){
  # calculate average abundance per species by site/year
  # adjusted by the largest average
  site_year_abund <- data %>% select(year, siteID, species) %>% 
    group_by(species, year, siteID) %>% 
    summarise(count = n())
  avg_abund <- site_year_abund %>% 
    group_by(species) %>% 
    summarise(avg_abund = mean(count)) %>% 
    mutate(adj_avg_abund = avg_abund/max(avg_abund))
  return(avg_abund)
}

find_reg_persist <- function(data){
  # create column for relative regional persistence (years overall)
  dat <- select(data, species, year) %>% 
    group_by(species, year) %>% 
    summarise(year_count = n())
  dat1 <- select(dat, species) %>% 
    group_by(species) %>% 
    summarise(years = n())
  total_years <- length(unique(data$year))
  dat1 <- mutate(dat1, rel_reg_persist = years/total_years)
  return(dat1)
}

find_reg_occup <- function(data){
  # create column for relative regional occupancy (sites overall)
  dat <- select(data, species, siteID) %>% 
    group_by(species, siteID) %>% 
    summarise(site_count = n())
  dat1 <- select(dat, species) %>% 
    group_by(species) %>% 
    summarise(sites = n())
  total_sites <- length(unique(data$siteID))
  dat1 <- mutate(dat1, rel_reg_occup = sites/total_sites)
  return(dat1)
}

join_regionals <- function(data){
  # combine the relative regional persistance and regional occupancy into one data frame
  dat <- left_join(x = find_reg_persist(data), y = find_reg_occup(data), by = "species")
  return(dat)
}

rel_local_occup <- function(data){
  # create a column for average relative local occupancy (site by year)
  dat <- select(data, year, species, siteID) %>% 
    group_by(year, species) %>% 
    summarise(site_by_year = n_distinct(siteID))
  total_sites <- length(unique(data$siteID))
  dat <- mutate(dat, rel_sites_by_year = site_by_year/total_sites)
  dat1 <- dat %>% group_by(species) %>% 
    summarise(mean_local_occup = mean(rel_sites_by_year))
}

rel_local_persist <- function(data){
  # create a column for average relative local persistance (year by site)
  dat <- select(data, year, species, siteID) %>% 
    group_by(siteID, species) %>% 
    summarise(years_by_site = n_distinct(year))
  total_years <- length(unique(data$year))
  dat <- mutate(dat, rel_years_by_site = years_by_site/total_years)
  dat1 <- dat %>% group_by(species) %>% 
    summarize(mean_local_persist = mean(rel_years_by_site))
  return(dat1)
}

join_locals <- function(data){
  # combine average relative local persistance and occupancy into one data frame
  dat <- left_join(x = rel_local_persist(data), y = rel_local_occup(data), by = "species")
  return(dat)
}

combine_all <- function(data){
  # combine regional and local persistance and occupancy into one data frame
  dat <- left_join(x = join_regionals(data), y = join_locals(data), by = "species")
  return(dat)
}

all_together <- function(data){
  # run all functions together to get the full output
  dat <- select_data(data) # can use pipes here
  dat <- group_data(dat)
  dat <- add_siteID(dat)
  dat <- combine_all(dat)
  return(dat)
}