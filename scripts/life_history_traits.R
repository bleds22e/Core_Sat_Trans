# Life History Traits - Rodents
  # for Multivariate stats project

#######################
# LIBRARIES

library(dplyr)

#######################
# DATA

pantheria <- read.table("data/Rodents/PanTHERIA.txt", header = TRUE, sep = "\t", na.strings = c("-999", "-999.00"))
pantheria <- dplyr::select(pantheria, MSW05_Order, MSW05_Genus, MSW05_Species, MSW05_Binomial, X5.1_AdultBodyMass_g, X6.1_DietBreadth, X12.1_HabitatBreadth, X22.1_HomeRange_km2) %>% 
             dplyr::arrange(MSW05_Genus)
pantheria <- rename(pantheria, Order = MSW05_Order, Genus = MSW05_Genus, Species = MSW05_Species, Binomial = MSW05_Binomial, 
                    AdultBodyMass_g = X5.1_AdultBodyMass_g, DietBreadth = X6.1_DietBreadth, HabitatBreadth = X12.1_HabitatBreadth,
                    HomeRange_km2 = X22.1_HomeRange_km2)
write.csv(pantheria, "data/Rodents/pantheria_all.csv")
write.csv(panth)
write.csv(all_data, "data/Rodents/all_data.csv")
