
localLinearSummaryGP <- function(subDf, Ind, shape, Nnew, x, y,
                                 sigma2Samples, terms, Name = NULL) {

  require(mvtnorm)

  ## Components of data matrix
  shapeNewPoints <- newPoints(Nnew, subDf, shape)
  dataSub <- x[Ind, ]
  dataSubx <- dataSub[, 1:3]

  meanSub <- colMeans(dataSubx)
  covSub <- cov(dataSubx)

  ## Create new data
  subNewPoints <- newPoints(Nnew, subDf, shape) #Points

  newx <- rmvnorm(Nnew, meanSub, covSub)
  
  newxsp <- newx %>% cbind(subNewPoints) %>% as.matrix()

  prednew <- predictGPposterior(newxsp, x, y, sigma2Samples, GPfit$hyper, Ke)

  ## Posterior projections
  xp <- as.matrix(cbind(1, newx))
  XtXinv <- solve(crossprod(xp))
  coefPost <- XtXinv %*% crossprod(xp, prednew)

  ## Credible intervals
  termsInt <- c("(Intercept)", terms)
  
  ci <- posteriorCI(coefPost, termsInt) %>%
    mutate(Name = Name)

  gamma <- xp %*% coefPost

  newrsq <- rsqGamma(gamma, prednew)

  return(list(
    meanSub = meanSub,
    covSub = covSub,
    newPoints = subNewPoints,
    newx = newx,
    ci = ci,
    gamma = gamma,
    rsqGamma = data.frame(rsq = newrsq, Name = Name)
  ))

}
