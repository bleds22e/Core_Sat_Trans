# Attempting to find all possible matrices
# 7/7/2016
# EKB

install.packages("permute")

m <- 9
n <- 5
hja_matrix <- matrix(sample(0:1, m * n, replace = TRUE), m, n)


