# Finding range of values per matrix
# 7/12/2016
# EKB, with math help from D. Harris

library(dplyr)
library(stringr)

yt_sev = 20
yt_jor = 13
yt_sgs = 7
yt_hja = 5

st_sev = 9
st_jor = 6
st_sgs = 7
st_hja = 9

find_range <- function(matrix){
  num <- sum(rowMeans(matrix)) / sum(rowSums(matrix) > 0)
  return(num)
}

#sum(rowMeans(matrix)) / sum(rowSums(matrix) > 0))
vec_max <- c(1,1,1,1,1,1,1)
vec_min <- c(0,0,0,0,0,0,1)

sgs_yp <- matrix(0, nrow = 7, ncol = 7)
sgs_yp1_max <- sgs_yp %>% replace(sgs_yp[1,], values = vec_max)
sgs_yp[1,] <- vec_max 
sgs_yp1_max <- sgs_yp
sgs_yp[2,] <- vec_max
sgs_yp2_max <- sgs_yp
sgs_yp[3,] <- vec_max
sgs_yp3_max <- sgs_yp
sgs_yp[4,] <- vec_max
sgs_yp4_max <- sgs_yp
sgs_yp[5,] <- vec_max
sgs_yp5_max <- sgs_yp
sgs_yp[6,] <- vec_max
sgs_yp6_max <- sgs_yp
sgs_yp[7,] <- vec_max
sgs_yp7_max <- sgs_yp

sgs_yp <- matrix(0, nrow = 7, ncol = 7)
sgs_yp1_min <- sgs_yp %>% replace(sgs_yp[1,], values = vec_min)
sgs_yp[1,] <- vec_min 
sgs_yp1_min <- sgs_yp
sgs_yp[2,] <- vec_min
sgs_yp2_min <- sgs_yp
sgs_yp[3,] <- vec_min
sgs_yp3_min <- sgs_yp
sgs_yp[4,] <- vec_min
sgs_yp4_min <- sgs_yp
sgs_yp[5,] <- vec_min
sgs_yp5_min <- sgs_yp
sgs_yp[6,] <- vec_min
sgs_yp6_min <- sgs_yp
sgs_yp[7,] <- vec_min
sgs_yp7_min <- sgs_yp

find_range(sgs_yp1_max)
find_range(sgs_yp2_max)
find_range(sgs_yp3_max)
find_range(sgs_yp4_max)
find_range(sgs_yp5_max)
find_range(sgs_yp6_max)
find_range(sgs_yp7_max)
find_range(sgs_yp1_min)
find_range(sgs_yp2_min)
find_range(sgs_yp3_min)
find_range(sgs_yp4_min)
find_range(sgs_yp5_min)
find_range(sgs_yp6_min)
find_range(sgs_yp7_min)

#### HJA

vec_max <- c(1,1,1,1,1,1,1,1,1)
vec_min <- c(0,0,0,0,0,0,0,0,1)

hja_yp <- matrix(0, nrow = 5 , ncol = 9)
hja_yp1_max <- hja_yp %>% replace(hja_yp[1,], values = vec_max)
hja_yp[1,] <- vec_max 
hja_yp1_max <- hja_yp
hja_yp[2,] <- vec_max
hja_yp2_max <- hja_yp
hja_yp[3,] <- vec_max
hja_yp3_max <- hja_yp
hja_yp[4,] <- vec_max
hja_yp4_max <- hja_yp
hja_yp[5,] <- vec_max
hja_yp5_max <- hja_yp

hja_yp <- matrix(0, nrow = 5 , ncol = 9)
hja_yp1_min <- hja_yp %>% replace(hja_yp[1,], values = vec_min)
hja_yp[1,] <- vec_min 
hja_yp1_min <- hja_yp
hja_yp[2,] <- vec_min
hja_yp2_min <- hja_yp
hja_yp[3,] <- vec_min
hja_yp3_min <- hja_yp
hja_yp[4,] <- vec_min
hja_yp4_min <- hja_yp
hja_yp[5,] <- vec_min
hja_yp5_min <- hja_yp


################################################
