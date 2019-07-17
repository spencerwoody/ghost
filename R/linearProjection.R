

linearProjection <- function(X, betaSamples, selectedVars) {

  Xbeta <- X %*% betaSamples
  
  Xred <- X[, selectedVars]
  
  betaDSS <- solve(crossprod(Xred), crossprod(Xred, Xbeta))

  betaDSSall <- matrix(nrow = nrow(betaSamples), ncol = ncol(betaSamples))

  betaDSSall[selectedVars, ] <- betaDSS
  betaDSSall[-selectedVars, ] <- 0

  return(betaDSSall)
  
}
