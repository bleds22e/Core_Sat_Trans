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

r <- 5
c <- 6
m1 <- round(matrix(runif(r*c), r, c))
m1
sum(rowSums(m1))
sum(colSums(m1))

matdensity <- 45
matsize <- 10
mat45 <- matrix(sample(c(rep(1, matdensity), rep(0, matsize*2 -
                                                   matdensity))), matsize, matsize)


mat <- matrix(0, 5, 6)
tofill <- sample(1:100)
for(i in 1:length(tofill)) {
  mat[tofill[i]] <- 1
}

tofill
