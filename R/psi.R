
psi <- function(X, betaSamples, betaProj, sigma2Samples) {

  if (is.vector(betaProj)) {
    betaProj <- matrix(rep(betaProj, ncol(betaSamples)), ncol = ncol(betaSamples))
  }

  yhatSamples <- X %*% betaSamples
  gammaSamples <- X %*% betaProj

  sqrt(
    colMeans((yhatSamples - gammaSamples)^2) + sigma2Samples
  ) - sqrt(sigma2Samples)       

}


