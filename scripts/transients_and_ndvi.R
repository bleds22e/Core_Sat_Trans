# Portal Transients and Productivity
# August 2018
# EKB

# LIBRARIES #
library(tidyverse)
library(portalr)

# DATA #

rodents <- read.csv("~bleds22e/Dropbox (UFL)/Git/PortalData/Rodents/Portal_rodent.csv", 
                    na.strings = "", stringsAsFactors = FALSE)
ndvi <- read.csv('~bleds22e/Dropbox (UFL)/Git/PortalData/NDVI/monthly_NDVI.csv', stringsAsFactors = FALSE)
