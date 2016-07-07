# make a graph of a hypothetical clustered distribution
  # with bimodal histograms on the side

hypothetical <- read.csv("hypotheticals.csv")

install.packages("ggExtra")
library(ggplot2)

plot <- ggplot(hypothetical, aes(x = Cluster_Occ, y = Cluster_Pers)) +
              geom_point(size = 3) +
              xlab("Regional Occupancy") +
              ylab("Regional Persistence") +
              theme_classic()
plot
ggExtra::ggMarginal(plot, type = "histogram", binwidth = 0.1)

dev.copy(png, "clustered_hist.png")
dev.off()

ggsave(file = "clustered_with_histograms.png", width = 6, height = 5)

ggsa