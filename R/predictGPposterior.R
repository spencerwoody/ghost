
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param datanew 
##' @param x original data
##' @param y vector of responses
##' @param sigma2Samples 
##' @param hyper Estimated hyperparameters from GP
##' @param Keigen Eigendecomposition of covariance matrix of f
##' @param verbose if TRUE (default) print out progress message
##' @return 
##' @author Spencer Woody
predictGPposterior <- function(datanew, x, y, sigma2Samples, hyper,
                               Keigen, verbose = TRUE) {
  

  require(GPFDA)
  require(dplyr)
  ## require(doParallel)

  ## registerDoParallel(makeCluster(detectCores()))

  N <- length(y)
  Nnew <- nrow(datanew)
  NMC <- length(sigma2Samples)

  hyperest <- hyper
  hyperest$vv <- NULL

  sigma2hat <- mean(sigma2Samples)

  ## Create object for output
  predMat <- matrix(nrow = Nnew, ncol = NMC)

  Kjoint <- cov.pow.ex(hyperest, Data = x, Data.new = datanew) +
    cov.linear(hyperest, Data = x, Data.new = datanew)

  Kstarstar <- cov.pow.ex(hyperest, Data = datanew) +
    cov.linear(hyperest, Data = datanew)

  KV <- Kjoint %*% Keigen$vectors
  Vty <- crossprod(Keigen$vectors, y)

  KcondEigen <- 1/(Keigen$values + sigma2hat)

  yhatmat <- KV %*% Diagonal(N, KcondEigen)
  yhatnew <- as.numeric(yhatmat %*% Vty)

  yhatcov <- as.matrix(Kstarstar -
                       yhatmat %*% t(Keigen$vectors) %*% t(Kjoint))

  predMat <- rmvnorm(NMC, yhatnew, yhatcov)

  ## output
  t(predMat) 

}

