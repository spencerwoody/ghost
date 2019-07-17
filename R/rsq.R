
##' R-squared posterior
##'
##' Calculates R-squared for posterior fit
##' @title 
##' @param y n-vector response
##' @param fitmat n \times NMC fitted values; or n-vector of single fit
##' @return 
##' @author Spencer Woody

rsq <- function(y, fitmat) {

  if (!is.matrix(fitmat)) {
    fitmat <- matrix(fitmat, ncol = 1)
  }

  SST <- sum((y - mean(y))^2)
  
  out <- rep(NA, ncol(fitmat))

  for (j in 1:ncol(fitmat)) {
    out[j] <- 1 - sum((y - fitmat[, j])^2) / SST
  }

  return(out)
}


