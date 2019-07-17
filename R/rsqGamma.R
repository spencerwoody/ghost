
rsqGamma <- function(gamma, yhat) {

  if (!is.matrix(gamma)) {
    gamma <- matrix(gamma, ncol = 1)
  }

  if (!is.matrix(yhat)) {
    yhat <- matrix(yhat, ncol = 1)
  }

  yhatBar <- colMeans(yhat)

  SST <- colSums(sweep(yhat, 2, yhatBar)^2)

  if (ncol(gamma) == 1) {
    gammaMat <- matrix(rep(as.numeric(gamma), ncol(yhat)),
                       ncol = ncol(yhat))

    SSR <- colSums((yhat - gammaMat)^2)
  } else {
    SSR <- colSums((yhat - gamma)^2)
  }

  1 - SSR / SST

}
