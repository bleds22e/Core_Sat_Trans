# Core-satellite-transient with Portal data
# Feb 17, 2016
# Ellen Bledsoe

# download libraries
library(dplyr)
library(ggplot2)

# read files

species <- read.csv("data/species_portal.csv", header = TRUE, na.strings = c(""))
surveys <- read.csv("data/surveys.csv", header = TRUE, na.strings = c(""))

# get only known rodents
rodents <- select(species, New.Code, Rodent) %>% 
           filter(Rodent == 1)
known_rodents <- select(rodents, New.Code) %>% 
                 filter(New.Code != 'DX', New.Code != 'NX',
                        New.Code != 'OX', New.Code != 'PX',
                        New.Code != 'RX', New.Code != 'SX', New.Code != 'UR')
colnames(known_rodents)[1] <- "species"

# get relevant survey periods & only control plots
surveys_control <- surveys %>% 
                   filter(period > 0 & period <= 436,
                          plot == '2' | plot == '4' | plot == '8' | plot == '11' |
                          plot == '12' | plot == '14' | plot == '17' | plot == '22')

portal <- semi_join(surveys_control, known_rodents, by = "species")

### calculating relative years and relative sites 

colnames(portal)[4] <- "year"
port_data <- select(portal, year, plot, species) %>% 
             group_by(species, year, plot) %>% 
             arrange(species, year, plot) %>% 
             summarise(count = n())

# number of years in which species occur overall

port_years <- select(port_data, species, year) %>% 
              group_by(species, year) %>% 
              summarise(year_count = n())
port_years_count <- select(port_years, species) %>% 
                    group_by(species) %>% 
                    summarise(years = n())

total_years <- length(unique(port_data$year))

port_years_count <- mutate(port_years_count, rel_years = years/total_years)

# number of sites in which species occur overall

port_plots <- select(port_data, species, plot) %>% 
              group_by(species, plot) %>% 
              summarise(plot_count = n())
port_plots_count <- select(port_plots, species) %>% 
                    group_by(species) %>% 
                    summarise(plots = n())

total_plots <- length(unique(port_data$plot))

port_plots_count <- mutate(port_plots_count, rel_sites = plots/total_plots)

# average relative number of plots per year

port_data1 <- select(portal, year, plot, species) %>% 
              group_by(species, year, plot) %>% 
              arrange(species, year, plot) 

port_plots_by_year <- select(port_data1, year, species, plot) %>% 
                      group_by(year, species) %>% 
                      summarise(plot_by_year = n_distinct(plot))

port_rel_plots_by_year <- mutate(port_plots_by_year, rel_plots_by_year = plot_by_year/total_plots)

port_plots_year_species <- port_rel_plots_by_year %>% 
                           group_by(species) %>% 
                           summarize(mean_rel_sites = mean(rel_plots_by_year))

# average relative number of years per plot

port_years_by_plot <- select(port_data1, year, species, plot) %>% 
                      group_by(plot, species) %>% 
                      summarise(years_by_plot = n_distinct(year))

port_rel_years_by_plot <- mutate(port_years_by_plot, rel_years_by_plot = years_by_plot/total_years)

port_years_plot_species <- port_rel_years_by_plot %>% 
                           group_by(species) %>% 
                           summarize(mean_rel_years = mean(rel_years_by_plot))

### combine all relevant data into one dataframe

portal_rel_data <- inner_join(port_years_count, port_plots_count, by = "species")
portal_mean_rel_data <- inner_join(port_years_plot_species, port_plots_year_species, by = "species")
portal_total_rel_data <- inner_join(portal_rel_data, portal_mean_rel_data, by = "species")

### make scatterplots of portal control data

ggplot(portal_total_rel_data, aes(x = rel_sites, y = rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied") + 
  ylab("Relative Number of Years Present") +
  ggtitle("Portal: Control Only") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal_broad_yr_site.png", width = 6, height = 5)

ggplot(portal_total_rel_data, aes(x = mean_rel_sites, y = rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied Per Year") + 
  ylab("Relative Number of Years Present") +
  ggtitle("Portal: Control Only") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal_mean_rel_site.png", width = 6, height = 5)

ggplot(portal_total_rel_data, aes(x = rel_sites, y = mean_rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied") + 
  ylab("Relative Number of Years Present Per Site") +
  ggtitle("Portal: Control Only") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal_mean_rel_year.png", width = 6, height = 5)

ggplot(portal_total_rel_data, aes(x = mean_rel_sites, y = mean_rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied Per Site") + 
  ylab("Relative Number of Years Present Per Site") +
  ggtitle("Portal: Control Only") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal_mean_rel_year&site.png", width = 6, height = 5)

### everything all over except with long-term krat exclosures included

# get relevant survey periods & only control plots
surveys_control_krat <- surveys %>% 
                        filter(period > 0 & period <= 436,
                               plot == '2' | plot == '4' | plot == '8' | plot == '11' |
                               plot == '12' | plot == '14' | plot == '17' | plot == '22' |
                               plot == '3' | plot == '15' | plot == '19' | plot == '21')

portal1 <- semi_join(surveys_control_krat, known_rodents, by = "species")

# calculating relative years and relative sites 

colnames(portal1)[4] <- "year"
port_data1 <- select(portal1, year, plot, species) %>% 
              group_by(species, year, plot) %>% 
              arrange(species, year, plot) %>% 
              summarise(count = n())

# number of years in which species occur overall

port_years1 <- select(port_data1, species, year) %>% 
               group_by(species, year) %>% 
               summarise(year_count = n())
port_years_count1 <- select(port_years1, species) %>% 
                     group_by(species) %>% 
                     summarise(years = n())

total_years <- length(unique(port_data$year))

port_years_count1 <- mutate(port_years_count1, rel_years = years/total_years)

# number of sites in which species occur overall

port_plots1 <- select(port_data1, species, plot) %>% 
               group_by(species, plot) %>% 
               summarise(plot_count = n())
port_plots_count1 <- select(port_plots1, species) %>% 
                     group_by(species) %>% 
                     summarise(plots = n())

total_plots1 <- length(unique(port_data1$plot))

port_plots_count1 <- mutate(port_plots_count1, rel_sites = plots/total_plots1)

# average relative number of plots per year

port_data1 <- select(portal1, year, plot, species) %>% 
              group_by(species, year, plot) %>% 
              arrange(species, year, plot) 

port_plots_by_year1 <- select(port_data1, year, species, plot) %>% 
                       group_by(year, species) %>% 
                       summarise(plot_by_year = n_distinct(plot))

port_rel_plots_by_year1 <- mutate(port_plots_by_year1, rel_plots_by_year = plot_by_year/total_plots1)

port_plots_year_species1 <- port_rel_plots_by_year1 %>% 
                            group_by(species) %>% 
                            summarize(mean_rel_sites = mean(rel_plots_by_year))

# average relative number of years per plot

port_years_by_plot1 <- select(port_data1, year, species, plot) %>% 
                       group_by(plot, species) %>% 
                       summarise(years_by_plot = n_distinct(year))

port_rel_years_by_plot1 <- mutate(port_years_by_plot1, rel_years_by_plot = years_by_plot/total_years)

port_years_plot_species1 <- port_rel_years_by_plot1 %>% 
                            group_by(species) %>% 
                            summarize(mean_rel_years = mean(rel_years_by_plot))

### combine all relevant data into one dataframe

portal_rel_data1 <- inner_join(port_years_count1, port_plots_count1, by = "species")
portal_mean_rel_data1 <- inner_join(port_years_plot_species1, port_plots_year_species1, by = "species")
portal_total_rel_data1 <- inner_join(portal_rel_data1, portal_mean_rel_data1, by = "species")

### make scatterplots of portal control data

ggplot(portal_total_rel_data1, aes(x = rel_sites, y = rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied") + 
  ylab("Relative Number of Years Present") +
  ggtitle("Portal: Control & Krat Exclosure") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal1_broad_yr_site.png", width = 6, height = 5)

ggplot(portal_total_rel_data1, aes(x = mean_rel_sites, y = rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied Per Year") + 
  ylab("Relative Number of Years Present") +
  ggtitle("Portal: Control & Krat Exclosure") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal1_mean_rel_site.png", width = 6, height = 5)

ggplot(portal_total_rel_data1, aes(x = rel_sites, y = mean_rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied") + 
  ylab("Relative Number of Years Present Per Site") +
  ggtitle("Portal: Control & Krat Exclosure") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal1_mean_rel_year.png", width = 6, height = 5)

ggplot(portal_total_rel_data1, aes(x = mean_rel_sites, y = mean_rel_years)) +
  geom_point(size = 3) +
  xlab("Relative Number of Sites Occupied Per Site") + 
  ylab("Relative Number of Years Present Per Site") +
  ggtitle("Portal: Control & Krat Exclosure") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_line(colour = "light gray"))
ggsave(file = "portal1_mean_rel_year&site.png", width = 6, height = 5)
