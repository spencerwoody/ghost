
library(here)
library(ggplot2)
library(dplyr)
library(rpart)
library(rpart.plot)
library(latex2exp)
library(GPFDA)
library(Matrix)
library(mvtnorm)

load("R-output/CA-GPfit.Rdata") 

Rcpp::sourceCpp(here("src/multiplyColumns.cpp"))

###############################################################################
                                        #         Precache computation        #
###############################################################################

NMC <- 1000

yhat <- GPfit$fitted.mean

KeVec <- Ke$vectors
KeVal <- Ke$values

KeVect <- t(Ke$vectors)

sigma2hat <- exp(GPfit$hyper$vv)

KpostEigen <- 1 / (1 / Ke$values + 1 / sigma2hat)
KpostHalf <- multiplyColumns(Ke$vectors, sqrt(KpostEigen))

## Kpost <- tcrossprod(KpostHalf, KpostHalf)

yhat21 <- crossprod(KpostHalf, y)

KpostHalf <- KeVec %*% Diagonal(N, KpostEigen)

yhat21 <- t(KpostHalf) %*% y

yhat2 <- KeVec %*% yhat21 / sigma2hat

cor(yhat, yhat2 %>% as.numeric())

###############################################################################
                                        #       Posterior samples of GP       #
###############################################################################

nburn <- 100

yhatSamples <- matrix(nrow = N, ncol = NMC + nburn + 1)
sigma2Samples <- rep(NA, NMC + nburn + 1)

yhatSamples[, 1] <- yhat
sigma2Samples[1] <- sigma2hat

prog <- progress_estimated(NMC + nburn)

for (j in 2:(NMC + nburn + 1)) {

  RSSj <- sum((y - yhatSamples[, j - 1])^2)

  sigma2Samples[j] <- 1 / rgamma(1, N / 2, RSSj / 2)

  ## Covariance matrix squareroot
  KpostEigenj <- 1 / (1 / Ke$values + 1 / sigma2Samples[j])
  KpostHalfj <- multiplyColumns(Ke$vectors, sqrt(KpostEigenj))

  ## Posterior mean
  yhatj <- KpostHalfj %*% crossprod(KpostHalfj, y / sigma2Samples[j])
  
  ## Draw from conditional posterior
  yhatSamples[, j] <- KpostHalfj %*% rnorm(N) + yhatj
  
  ## if (j %% 100 == 0) cat(sprintf("%i out of %i...\n", j, NMC + nburn))
  prog$tick()$print()

  cat("iteration %i out of %i...\n", j, NMC + nburn)
  
}

yhatSamples <- yhatSamples[, -(1:(nburn+1))]
sigma2Samples <- sigma2Samples[-(1:(nburn+1))]

save(yhatSamples, sigma2Samples, file = "R-output/GPgibbs-samples.Rdata")

