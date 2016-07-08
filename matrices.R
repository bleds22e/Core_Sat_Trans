# Attempting to find all possible matrices
# 7/7/2016
# EKB

install.packages("permute")
library(permute)

m <- 6
n <- 5
hja_matrix <- matrix(sample(0:1, m * n, replace = TRUE), m, n)
hja_matrix
rowSums(hja_matrix)
colSums(hja_matrix)
