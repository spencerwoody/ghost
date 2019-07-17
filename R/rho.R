##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param X N x p
##' @param y 
##' @param betaSamples p x NMC
##' @param gammaSamples 
##' @param sigma2Samples 
##' @return 
##' @author Spencer Woody
rho <- function(X, betaSamples, betaProj, sigma2Samples) {

  if (is.vector(betaSamples)) {
    betaSamples <- matrix(rep(betaProj, length(sigma2Samples)), ncol = length(sigma2Samples))
  }

  if (is.vector(betaProj)) {
    betaProj <- matrix(rep(betaProj, ncol(betaSamples)), ncol = ncol(betaSamples))
  }

  yhatSamples <- X %*% betaSamples
  gammaSamples <- X %*% betaProj
  
  numerator <- colMeans(yhatSamples^2)

  denominator <- numerator + sigma2Samples +
    colMeans((yhatSamples - gammaSamples)^2)

  numerator / denominator

}
