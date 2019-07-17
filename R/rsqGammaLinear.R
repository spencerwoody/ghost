
rsqGammaLinear <- function(X, betaSamples, betaProj) {
  yhat <- X %*% betaSamples
  gamma <- X %*% betaProj

  rsqGamma(gamma, yhat)

} 
