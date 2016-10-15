setwd('~/Github/BAC/')
source('CalcLogLike_function.R')
source('LogPriorOdds_function.R')
source('UpdateAlphaX_function.R')
source('UpdateAlphaY_function.R')
source('UpdateBeta_function.R')
source('BAC_function.R')

num_conf <- 20
N <- 300
D <- rnorm(N * num_conf, mean = 0, sd = 1)
D <- matrix(D, nrow = N, ncol = num_conf)
X <- D[, 1] + D[, 2] + D[, 3] + rnorm(N, 0, 1)
Y <- D[, 1] + D[, 2] + rnorm(N, 0, 1)

bac <- BAC(X = X, Y = Y, D = D, Nsims = 1000)

apply(bac$alphas$X, 2, mean)
apply(bac$alphas$Y, 2, mean)
mean(bac$beta)
lm(Y ~ D[, 1] + D[, 2])$coef[2]
