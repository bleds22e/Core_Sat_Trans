# Figure for Transients Framework Paper
# EKB
# July-August 2020

# LIBRARIES #

library(ggplot2)
library(patchwork)

# PLOTS #

x <- as.data.frame(cos(seq(1*pi, 3*pi, 0.01)))
y <- as.data.frame(cos(seq(0.75*pi, 2.75*pi, 0.01)))
n <- as.data.frame(1:nrow(x))
data <- bind_cols(x, y, n)
colnames(data) <- c("Carrying Capacity",
                    "Population",
                    "Index")

data <- pivot_longer(data, 1:2, names_to = "var")


(plot1 <- ggplot(filter(data, var == "Carrying Capacity"), aes(Index, value)) + 
    geom_hline(yintercept = 0, color = "gray") +
    geom_line(aes(linetype = var)) +
    scale_y_continuous(limits = c(-1,1), breaks = c(0), 
                       labels = c("Avg."), sec.axis = dup_axis()) +
    ylab("NDVI") +
    theme_classic() +
    theme(axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1),
          legend.position = "none"))

(plot2 <- ggplot(data, aes(Index, value)) + 
    geom_hline(yintercept = 0, color = "gray") +
    geom_line(aes(linetype = var)) +
    scale_y_continuous(limits = c(-1,1), breaks = c(0), 
                       labels = c("Avg."), sec.axis = dup_axis()) +
    ylab("Number of Individuals") +
    xlab("Time") +
    theme_classic() +
    theme(axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1),
          legend.title = element_blank()))

(multipanel <- (plot1/plot2))

ggsave("plots/fig1.png")
