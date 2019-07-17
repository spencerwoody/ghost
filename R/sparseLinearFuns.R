


linearProjection <- function(X, betaSamples, selectedVars) {

  Xbeta <- X %*% betaSamples
  
  Xred <- X[, selectedVars]
  
  betaDSS <- solve(crossprod(Xred), crossprod(Xred, Xbeta))

  betaDSSall <- matrix(nrow = nrow(betaSamples), ncol = ncol(betaSamples))

  betaDSSall[selectedVars, ] <- betaDSS
  betaDSSall[-selectedVars, ] <- 0

  return(betaDSSall)
  
}


getLargestBetas <- function(beta, lambda) {
  ## Using beta and lambda outputs from glmnet to get point estimates of DSS

  require(dplyr)
  
  betaLambda <- NULL
  lambdaVec <- NULL
  modelList <- list()

  count <- 1

  for (j in rev(1:ncol(beta))) {

    nonzeroIdx <- (beta[, j] != 0)

    lengthI <- sum(!nonzeroIdx)

    if (lengthI == ncol(X)) {
      next
    }
    
    if (j == ncol(beta)) {

      betaLambda <- cbind(beta[, j], betaLambda)
      lambdaVec <- c(lambda[j], lambdaVec)
      
      modelList[[count]] <- which(nonzeroIdx) %>% as.numeric()

      count <- count + 1
      
      lengthOld <- lengthI
      
    } else if (lengthI != lengthOld) {

      betaLambda <- cbind(beta[, j], betaLambda)
      lambdaVec <- c(lambda[j], lambdaVec)
      lengthOld <- lengthI
      
      modelList[[count]] <- which(nonzeroIdx) %>% as.numeric()
      
      count <- count + 1

      
    } else {
      next
    }
    
  }

  
  return(list(
    betaLambda = betaLambda,
    lambda = lambdaVec,
    models = modelList
  ))
  
}

projectionTreatment <- function(Y, X, D, modelList, alphaSamples, betaSamples,
                                confLevel = 0.05, verbose = F) {

  nModels <- length(modelList)

  ## Dataframe for storing the projected models
  alphaDf <- NULL

  ## Dataframe for storing projected model summary
  alphaSummaryDf <- NULL

  ## List for storing entire projected posterior sample matix
  fullProjList <- vector("list", nModels)

  ## Posterior fitted values
  yHat <- X %*% betaSamples + alphaSamples * D
  
  for (j in 1:nModels) {
    
    ## Design matrix for this model
    Xd <- cbind(D, X[, modelList[[j]]])
    
    ## Projected samples
    fullProjList[[j]] <- solve(crossprod(Xd), crossprod(Xd, yHat))
    
    ## Obtain treatment effect estimate
    alphaProjSamples <- fullProjList[[j]][1, ] 
    
    ## Store these samples
    alphaDf <- alphaDf %>%
    rbind(data.frame(
      model = length(modelList[[j]]),
      alphaSamples = alphaProjSamples
    ))
    
    ## Calculate mean and 95% credible intervals
    alphaSummaryDf <- alphaSummaryDf %>%
    rbind(data.frame(
      postMean = mean(alphaProjSamples),
      CIlo = quantile(alphaProjSamples, 0.025),
      CIhi = quantile(alphaProjSamples, 0.975),
      model = length(modelList[[j]])
    ))

    if (verbose) print(j)

  }

  return(list(
    alphaDf = alphaDf,
    alphaSummaryDf = alphaSummaryDf,
    fullProjList = fullProjList
  ))
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param X 
##' @param betaSamples 
##' @param sigma2Samples 
##' @param betaLambda 
##' @return 
##' @author Spencer Woody
varExpFixed <- function(X, betaSamples, sigma2Samples, betaLambda) {
  Xbeta <- X %*% betaSamples

  XbetaNorm <- colMeans(Xbeta^2)

  XbetaLambda <- X %*% betaLambda

  excessError <- colMeans(
    sweep(Xbeta, 1, XbetaLambda, "-")^2
  )

  rho2 <- XbetaNorm / (XbetaNorm + sigma2Samples + excessError)

  return(rho2)
}

varExpPost <- function(X, betaSamples, sigma2Samples, betaDSS) {
  Xbeta <- X %*% betaSamples

  XbetaNorm <- colMeans(Xbeta^2)

  XbetaDSS <- X %*% betaDSS

  excessError <- colMeans((Xbeta - XbetaDSS)^2)
  
  rho2 <- XbetaNorm / (XbetaNorm + sigma2Samples + excessError)

  return(rho2)
}


getbetaDSS <- function(X, betaSamples, ix) {

  Xbeta <- X %*% betaSamples
  
  Xred <- X[, ix]
  
  betaDSS <- solve(crossprod(Xred), crossprod(Xred, Xbeta))

  betaDSSall <- matrix(nrow = nrow(betaSamples), ncol = ncol(betaSamples))

  betaDSSall[ix, ] <- betaDSS
  betaDSSall[-ix, ] <- 0

  return(betaDSSall)
  
}

