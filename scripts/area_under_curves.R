# Finding the Area Under the Time-series Curve
#     - transients and NDVI
# EKB
# Oct 2017

# LIBRARIES #
library(devtools)
library(tidyverse)
source('scripts/transientNDVI_exploration.R')

# Try different packages

# timeseriesr
install_github("aaronbenz/timeseriesr")
# StreamMetabolism
install.packages('StreamMetabolism')

# Temporary fix #

# from msProcess package
library(ifultools)
zeroCross <- function(x, slope="positive")
{
  ifultools::checkVectorType(x,"numeric")
  ifultools::checkScalarType(slope,"character")
  slope <- match.arg(slope,c("positive","negative"))

  ipost  <- splus2R::ifelse1(slope == "negative", sort(which(c(x, 0) < 0 & c(0, x) > 0)),
                    sort(which(c(x, 0) > 0 & c(0, x) < 0)))
  offset <- apply(matrix(abs(x[c(ipost-1, ipost)]), nrow=2, byrow=TRUE), MARGIN=2, order)[1,] - 2
  ipost + offset
}

NDVI_peak_ts <- ts(NDVI_peak$NDVIpeak, start = c(1992, 3), end = c(2014, 11), frequency = 12)
plot(NDVI_peak_ts)
abline(h = 0)
abline(h = sd(NDVI_peak$NDVIpeak), col = 'red')
abline(h = -sd(NDVI_peak$NDVIpeak), col = 'red')
       
zcross_pos <- zeroCross(NDVI_peak_ts, slope = 'positive')
zcross_neg <- zeroCross(NDVI_peak_ts, slope = 'negative')

# make linear equation between two closest points to each zero-crossing

find_x <- function(x, y){
  # finds x where y equals zero
  m <- (y[2]-y[1]) / (x[2]-x[1])
  b <- y[1] - (m*x[1]) 
  xout = (-b)/m
  return(xout)
}
find_x(x = c(-2,2), y = c(-2,2))

# add these points to the time series in the appropriate places
# then use maybe `simp()` or whatever to calculate area underneath
par(mfrow = c(2,1))

plot(NDVI_peak_ts)
abline(h = median(NDVI_peak$NDVIpeak))
abline(h = quantile(NDVI_peak$NDVIpeak, .75), col = 'red')
abline(h = quantile(NDVI_peak$NDVIpeak, .25), col = 'red')

plot(density((NDVI_peak$NDVIpeak)))
abline(v = median(NDVI_peak$NDVIpeak))
abline(v = quantile(NDVI_peak$NDVIpeak, .75), col = 'red')
abline(v = quantile(NDVI_peak$NDVIpeak, .25), col = 'red')


