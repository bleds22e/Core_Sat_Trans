# Making Theoretical Figure for Paper
# November 2019
# EKB

# Make individual figures then piece together
# Figures to make:
#   1. Pr(success|leave) vs. Productivity
#   2. Pr(success|stay) vs. Change in Productivity
#   3. Productivity vs. Time
#   4. Change in Productivity vs. Time
#   5. Leave|Stay vs. Time

# Central Tenant: Pr(success|leave) > Pr(success|stay)

library(tidyverse)
library(ggpubr)

# Plot 1: 

(plot1 <- ggplot(data.frame(x = c(0, 1.5)), aes(x)) + 
  stat_function(fun = function(x) exp(x^2) + 1) +
  xlab('Productivity') +
  ylab('P(success | leave)') +
  scale_y_continuous(limits = c(0, 10), 
                     breaks = c(0, 10),
                     labels = c(0.0, 1.0)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()))
#ggsave("plots/theoretical_plots/plot1.png", plot1)

# Plot 2: 

x_axis_title <- expression(Delta*" Productivity")

(plot2 <- ggplot(data.frame(x = c(-1.5, 1.5)), aes(x)) + 
  stat_function(fun = function(x) x^3 + x + 1) +
  xlab(x_axis_title) +
  ylab('P(success | stay)') +
  scale_y_continuous(limits = c(-5, 6), 
                     breaks = c(-5, 6),
                     labels = c(0.0, 1.0)) +
  scale_x_continuous(breaks = c(0), labels = c(0)) +  
  theme_classic()) 
#ggsave("plots/theoretical_plots/plot2.png", plot2)

# Plot 3: Productivity vs. Time

(plot3 <- ggplot(data.frame(x = c(-5, 5)), aes(x)) + 
  stat_function(fun = dnorm, args = list(sd = 1.5)) +
  xlab('Time') +
  ylab('Productivity') +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()))
#ggsave("plots/theoretical_plots/plot3.png", plot3)
  
# Plot 4: Change in Productivity vs. Time

y_axis_title <- expression(Delta*" Productivity")

(plot4 <- ggplot(data.frame(x = c(0, 2*pi)), aes(x)) + 
  stat_function(fun = sin) +
  geom_hline(yintercept = 0, color = "gray", linetype = "longdash") +
  scale_y_continuous(limits = c(-1,1), breaks = c(0), labels = c(0)) +
  xlab('Time') +
  ylab(y_axis_title) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()))
#ggsave("plots/theoretical_plots/plot4.png", plot4)

# Piece together
(multipanel <- ggarrange(plot3, plot4, plot1, plot2, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2))
#ggsave("plots/theoretical_plots/multipanel_figure.png", multipanel)
