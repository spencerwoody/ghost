

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param X 
##' @param y 
##' @param betaSamples 
##' @param sigma2Samples 
##' @param adaptive if TRUE (default), use adaptive lasso, weighting
##'   by the postserior mean
##' @param varnames 
##' @param alpha 
##' @param ... other arguments, e.g., to glmnet
##' @return
##' @author Spencer Woody
sparseLinearProj <- function(X, y, betaSamples, sigma2Samples,
                       adaptive = TRUE, varnames = NA, alpha = 0.05,
                       ...) {

  require(reshape2)
  require(dplyr)
  require(tidyr)
  require(lars)
  require(stringr)

  p <- ncol(X)

  if (any(is.na(varnames))) {
    varnames <- str_c("X", 1:p)
  }

  varnamesDf <- data.frame(
    Var1 = seq_along(varnames),
    varname = varnames
  )

  ## Posterior mean of fitted value for y (X %*% betabar)
  yhatSamples <- X %*% betaSamples

  yFit <- rowMeans(yhatSamples)

  if (adaptive == TRUE) {
    betaBar <- rowMeans(betaSamples)
    weights <- 1 / abs(betaBar)
    
    Xweighted <- sweep(X, 2, weights, "/")

    dss <- lars(Xweighted, yFit, intercept = F, normalize = F, ...)
    dss$beta <- t(dss$beta)
    dss$beta <- sweep(dss$beta, 1, weights, "/") %>% as.matrix()

  } else {
    dss <- glmnet(X, yFit, intercept = F, standardize = F)

    dssBetaFull <- dss$beta %>% as.matrix()

    dssLambdaFull <- dss$lambda 
  }

  ## coefficients
  betaLambda <- dss$beta %>% as.matrix()

  ## Remove all-zeros column
  if (all(betaLambda[, 1] == 0)) {
    betaLambda <- betaLambda[, -1]
  }

  ## Remove all-zeros column
  if (all(betaLambda[, ncol(betaLambda)] != 0)) {
    betaLambda <- betaLambda[, -ncol(betaLambda)]
  }

  ## Complexity parameters
  dssLambdaFull <- dss$lambda

  ## Selected variables (removing intecept-only model)
  selectedVars <- apply(betaLambda, 2, function(x) which(x != 0))

  ## Vector of numbers of coefficients in each model
  modelSize <- sapply(selectedVars, length) %>% as.numeric()

  ## Get projected linear estimates, with summarization statistics
  betaProjList <- vector("list", length(selectedVars))

  betaProjDfList <- vector("list", length(selectedVars))

  rhoList <- vector("list", length(selectedVars))
  psiList <- vector("list", length(selectedVars))
  rsqGammaList <- vector("list", length(selectedVars))
  phiGammaList <- vector("list", length(selectedVars))
  
  summaryDfList <- vector("list", length(selectedVars))


  for (k in 1:length(selectedVars)) {
    betaProjList[[k]] <- linearProjection(X, betaSamples, selectedVars[[k]])

    betaProjDfList[[k]] <- betaProjList[[k]] %>%
      melt() %>%
      dplyr::select(-Var2) %>%
      left_join(varnamesDf, by = "Var1") %>%
      mutate(modelSize = modelSize[k]) %>%
      dplyr::select(-Var1)
    
    rhoList[[k]] <- rho(X, betaSamples, betaProjList[[k]], sigma2Samples)
    psiList[[k]] <- psi(X, betaSamples, betaProjList[[k]], sigma2Samples)
    rsqGammaList[[k]] <- rsqGammaLinear(X, betaSamples, betaProjList[[k]])
    phiGammaList[[k]] <- phiGammaLinear(X, y, betaProjList[[k]], sigma2Samples)

    summaryDfList[[k]] <- data.frame(
      modelSize = modelSize[k],
      rho = rhoList[[k]],
      psi = psiList[[k]],
      phi_gamma = phiGammaList[[k]],
      rsq_gamma = rsqGammaList[[k]]
    )

  }

  modelSizeFull <- sprintf("Full (%i)", p)

  betaSamplesDf <- betaSamples %>%
    melt() %>%
      dplyr::select(-Var2) %>%
      left_join(varnamesDf, by = "Var1") %>%
      mutate(modelSize = modelSizeFull) %>%
      dplyr::select(-Var1)
  

  betaProjDf <- plyr::rbind.fill(betaProjDfList)

  modelSizeLevels <- betaProjDf %>%
    pull(modelSize) %>%
    as.numeric() %>%
    unique() %>% 
    sort(decreasing = TRUE)

  betaProjDf <- betaProjDf %>%
    rbind(betaSamplesDf) %>%
    mutate(modelSize = factor(modelSize, levels = c(modelSizeFull, modelSizeLevels)))
  
  varnameRank <- betaProjDf %>%
    group_by(varname) %>%
    summarize(Rank = mean(value == 0)) %>%
    arrange(Rank) %>%
    pull(varname)

  betaProjDf <- betaProjDf %>%
    mutate(varname = factor(varname, levels = varnameRank))

  betaProjSummary <- betaProjDf %>%
    group_by(modelSize, varname) %>%
    summarize(mid = mean(value),
              lo = quantile(value, alpha / 2),
              hi = quantile(value, 1 - alpha / 2))

  summaryDf <- plyr::rbind.fill(summaryDfList) %>%
    mutate(modelSize = factor(modelSize, levels = rev(unique(modelSize))))

  summaryDfLong <- summaryDf %>%
    gather(stat, value, -modelSize) 

  summaryDfCI <- summaryDfLong %>% 
    group_by(modelSize, stat) %>%
    summarize(mid = mean(value),
              lo = quantile(value, alpha / 2),
              hi = quantile(value, 1 - alpha / 2))

  ## summarydfCI <- summaryDf %>%
  ##   group_by(modelSize) %>%
  ##   summarize(lo = quantile())

  ## Pick out largest coefficient estimates for selected models (point
  ## estimates)
  lambdaVec <- NULL

  jVec <- NULL
  dfListCount <- 1

  dfList <- list()
  
  summaryDfOldList <- vector("list", ncol(betaLambda))
  eeDfList <- vector("list", ncol(betaLambda))

  betaLambdaDfList <- vector("list", ncol(betaLambda))

  for (j in 1:ncol(betaLambda)) {
    summaryDfOldList[[j]] <- data.frame(
      modelSize = sum(betaLambda[, j] != 0),
      rho = rho(X, betaSamples, betaLambda[, j], sigma2Samples),
      psi = psi(X, betaSamples, betaLambda[, j], sigma2Samples),
      rsq_gamma = rsqGammaLinear(X, betaSamples, betaLambda[, j]),
      phi_gamma = phiGammaLinear(X, y, betaLambda[, j], sigma2Samples)
    )

    betaLambdaDfList[[j]] <- betaLambda[, j] %>%
      matrix(ncol = 1) %>%
      melt() %>%
      dplyr::select(-Var2) %>%
      left_join(varnamesDf, by = "Var1") %>%
      mutate(modelSize = modelSize[j]) %>%
      dplyr::select(-Var1)
    
    eeDfList[[j]] <- data.frame(
      psi = mean(psi(X, betaSamples, betaLambda[, j], sigma2Samples)),
      rho = mean(rho(X, betaSamples, betaLambda[, j], sigma2Samples)),
      terms = varnames,
      coef = betaLambda[, j]
    )

  }

  betaLambdaDf <- plyr::rbind.fill(betaLambdaDfList) %>%
    mutate(modelSize = as.character(modelSize))

  summaryDfOld <- plyr::rbind.fill(summaryDfOldList) %>%
    mutate(modelSize = factor(modelSize,
                              levels = sort(unique(modelSize), decreasing = TRUE)))


  summaryDfOldLong <- summaryDfOld %>%
    gather(stat, value, -modelSize) 

  summaryDfOldCI <- summaryDfOldLong %>% 
    group_by(modelSize, stat) %>%
    summarize(mid = mean(value),
              lo = quantile(value, alpha / 2),
              hi = quantile(value, 1 - alpha / 2))


  eeDf <- plyr::rbind.fill(eeDfList)

  ## Output
  list(
    selectedVars = selectedVars,
    betaLambda = betaLambda, #Point estimates
    betaLambdaDf = betaLambdaDf, 
    betaProjList = betaProjList,
    betaProjDf = betaProjDf,
    betaProjSummary = betaProjSummary,
    rhoList = rhoList,
    psiList = psiList,
    rsqGammaList = rsqGammaList,
    phiGammaList = phiGammaList,
    summaryDf = summaryDf,
    summaryDfLong = summaryDfLong,
    summaryDfCI = summaryDfCI,
    summaryDfOld = summaryDfOld,
    summaryDfOldLong = summaryDfOldLong,
    summaryDfOldCI = summaryDfOldCI,
    eeDf = eeDf,
    modelSize = modelSize
  )

}
