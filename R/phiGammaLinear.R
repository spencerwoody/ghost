
phiGammaLinear <- function(X, y, betaProj, sigma2Samples) {

  gamma <- X %*% betaProj

  phiGamma(y, gamma, sqrt(sigma2Samples))

} 
