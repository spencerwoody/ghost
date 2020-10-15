

GPposterior <- function(y, GPfit, NMC = 1000, nburn = 100, verbose = TRUE) {

  require(GPFDA)
  require(dplyr)

  N <- length(y)

  K1 <- GPFDA::cov.pow.ex(GPfit$hyper, GPfit$train.x)
  K2 <- GPFDA::cov.linear(GPfit$hyper, GPfit$train.x)

  K <- K1 + K2
  
  Ke <- eigen(K)

  yhat <- GPfit$fitted.mean

  KeVec <- Ke$vectors
  KeVal <- Ke$values

  sigma2hat <- exp(GPfit$hyper$vv)

  ## Gibbs sampler

  yhatSamples <- matrix(nrow = N, ncol = NMC + nburn + 1)
  sigma2Samples <- rep(NA, NMC + nburn + 1)

  yhatSamples[, 1] <- yhat
  sigma2Samples[1] <- sigma2hat

  

  prog <- progress_estimated(NMC + nburn)
  
  for (j in 2:(NMC + nburn + 1)) {

    RSSj <- sum((y - yhatSamples[, j - 1])^2)

    sigma2Samples[j] <- 1 / rgamma(1, N / 2, RSSj / 2)

    ## Covariance matrix squareroot
    ## KpostEigenj <- 1 / (1 / KeVal + 1 / sigma2Samples[j])
    ## KpostHalfj <- multiplyColumns(Ke$vectors, sqrt(KpostEigenj))

    ## Posterior mean
    ## yhatj <- KpostHalfj %*% crossprod(KpostHalfj, y / sigma2Samples[j])

    ## Posterior covariance
    invMat <- solve(K + sigma2Samples[j] * diag(N))

    Covj <- sigma2Samples[j] * K %*% invMat

    ## Posterior mean
    Ej <- K %*% invMat %*% y

    ## Draw from conditional posterior
    ## yhatSamples[, j] <- KpostHalfj %*% rnorm(N) + yhatj

    yhatSamples[, j] <- rmvnorm(1, Ej, sigma = Covj)
    
    ## if (j %% 100 == 0) cat(sprintf("%i out of %i...\n", j, NMC + nburn))
    ## if (verbose) prog$tick()$print()

    if (verbose & (j - 1) %% 50 == 0) {
      cat(sprintf("iteration %i out of %i...\n", j - 1, NMC + nburn))
    }
    
  }

  yhatSamples <- yhatSamples[, -(1:(nburn+1))]
  sigma2Samples <- sigma2Samples[-(1:(nburn+1))]

  list(yhatmat = yhatSamples, sigma2Samples = sigma2Samples)

}
