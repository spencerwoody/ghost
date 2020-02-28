
library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)
library(mgcv)
library(stringr)
library(mvtnorm)
library(pracma)
library(Matrix)
library(GPFDA)
library(gridExtra)
library(cowplot)
library(ggridges)
library(mvtnorm)
library(rpart)
library(rpart.plot)
theme_set(theme_cowplot())

library(multibart)
library(dbarts)

source(here::here("R/gamProjFun.R"))
source(here::here("R/rsq.R"))
source(here::here("R/phiGamma.R"))
source(here::here("R/rsqGamma.R"))
source(here::here("R/GPposterior.R"))

###############################################################################
                                        #            Set parameters           #
###############################################################################

rho <- 0.5
rho2 <- 0.5

N <- 400
p <- 6

## Define matrix for covariates
Sigma <- matrix(0, nrow = p, ncol = p)

diag(Sigma) <- 1

Sigma[1, 3] <- rho2
Sigma[1, 4:6] <- c(rho, rho^2, rho^3)

Sigma[2, 3] <- rho2
Sigma[2, 4:6] <- c(rho, rho^2, rho^3)

Sigma[3, 4:6] <- c(rho, rho^2, rho^3)

Sigma[4, 5:6] <- c(rho, rho^2)

Sigma[5, 6] <- rho

for (i in p:1) {
  for (j in 1:i) {
    Sigma[i, j] <- Sigma[j, i]
  }
}

sigma <- sqrt(0.5)

###############################################################################
                                        #             simulations             #
###############################################################################

nsims <- 1000

XList <- vector("list", nsims)
fList <- vector("list", nsims)
yList <- vector("list", nsims)
GPfitList <- vector("list", nsims)
postList <- vector("list", nsims)

## set.seed(115)

## set.seed(1)

set.seed(1234)

for (k in 1:nsims) {

  X <- rmvnorm(N, sigma = Sigma)

  XList[[k]] <- X

  fList[[k]] <- sigmoid(- 2 * X[, 1] * X[, 2]) + (X[, 3] / 3)^3

  yList[[k]] <- fList[[k]] + rnorm(N, sd = sigma)

  GPfitList[[k]] <- gpr(XList[[k]], yList[[k]], Cov = "pow.ex")

  postList[[k]] <- GPposterior(yList[[k]], GPfitList[[k]], verbose = FALSE)

  cat(sprintf("Simulation %i out of %i...\n", k, nsims))

}

## save.image("newsim2.Rdata")

###############################################################################
                                        #            Make summaries           #
###############################################################################

## load("newsim2.Rdata")

XdfList <- vector("list", nsims)

gamList <- vector("list", nsims)

treeList <- vector("list", nsims)

rsqDfList <- vector("list", nsims)

rsqDfSummaryList <- vector("list", nsims)

for (k in 1:nsims) {

  yhat <- GPfitList[[k]]$fitted.mean
  yhatmat <- postList[[k]]$yhatmat

  ## rsqGamma(yhat, yList[[k]])

  gamListK <- vector("list", 16)
  
  XdfList[[k]] <- data.frame(
    x1 = XList[[k]][, 1],
    x2 = XList[[k]][, 2],
    x3 = XList[[k]][, 3],
    x4 = XList[[k]][, 4],
    x5 = XList[[k]][, 5],
    x6 = XList[[k]][, 6]
  )

  ## Additive GAM
  gam1 <- gam(yhat ~ s(x1) + s(x2) + s(x3) + s(x4) + s(x5) + s(x6),
                 data = XdfList[[k]])

  ## Summarization residuals and tree
  gamres <- yhat - gam1$fitted.values
  treeList[[k]] <- rpart(gamres ~ x1 + x2 + x3 + x4 + x5 + x6,
                         control = rpart.control(maxdepth=7), data = XdfList[[k]])

  ## GAM with one interaction
  gam12 <- gam(yhat ~ s(x1, x2) + s(x3) + s(x4) + s(x5) + s(x6), data = XdfList[[k]])
  gam13 <- gam(yhat ~ s(x1, x3) + s(x2) + s(x4) + s(x5) + s(x6), data = XdfList[[k]])
  gam14 <- gam(yhat ~ s(x1, x4) + s(x2) + s(x3) + s(x5) + s(x6), data = XdfList[[k]])
  gam15 <- gam(yhat ~ s(x1, x5) + s(x2) + s(x3) + s(x4) + s(x6), data = XdfList[[k]])
  gam16 <- gam(yhat ~ s(x1, x6) + s(x2) + s(x3) + s(x4) + s(x5), data = XdfList[[k]])
  gam23 <- gam(yhat ~ s(x2, x3) + s(x1) + s(x4) + s(x5) + s(x6), data = XdfList[[k]])
  gam24 <- gam(yhat ~ s(x2, x4) + s(x1) + s(x3) + s(x5) + s(x6), data = XdfList[[k]])
  gam25 <- gam(yhat ~ s(x2, x5) + s(x1) + s(x3) + s(x4) + s(x6), data = XdfList[[k]])
  gam26 <- gam(yhat ~ s(x2, x6) + s(x1) + s(x3) + s(x4) + s(x5), data = XdfList[[k]])
  gam34 <- gam(yhat ~ s(x3, x4) + s(x1) + s(x2) + s(x5) + s(x6), data = XdfList[[k]])
  gam35 <- gam(yhat ~ s(x3, x5) + s(x1) + s(x2) + s(x4) + s(x6), data = XdfList[[k]])
  gam36 <- gam(yhat ~ s(x3, x6) + s(x1) + s(x2) + s(x4) + s(x5), data = XdfList[[k]])
  gam45 <- gam(yhat ~ s(x4, x5) + s(x1) + s(x2) + s(x3) + s(x6), data = XdfList[[k]])
  gam46 <- gam(yhat ~ s(x4, x6) + s(x1) + s(x2) + s(x3) + s(x5), data = XdfList[[k]])
  gam56 <- gam(yhat ~ s(x5, x6) + s(x1) + s(x2) + s(x3) + s(x4), data = XdfList[[k]])

  ## Full projections

  gamListK[[1]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam1)
  gamListK[[2]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam12)
  gamListK[[3]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam13)
  gamListK[[4]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam14)
  gamListK[[5]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam15)
  gamListK[[6]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam16)
  gamListK[[7]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam23)
  gamListK[[8]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam24)
  gamListK[[9]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam25)
  gamListK[[10]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam26)
  gamListK[[11]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam34)
  gamListK[[12]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam35)
  gamListK[[13]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam36)
  gamListK[[14]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam45)
  gamListK[[15]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam46)
  gamListK[[16]] <- gamProjectionFun(yhatmat, XdfList[[k]], gam56)

  gam12 <- gam(yhat ~ s(x1, x2) + s(x3) + s(x4) + s(x5) + s(x6), data = XdfList[[k]])
  gam13 <- gam(yhat ~ s(x1, x3) + s(x2) + s(x4) + s(x5) + s(x6), data = XdfList[[k]])
  gam14 <- gam(yhat ~ s(x1, x4) + s(x2) + s(x3) + s(x5) + s(x6), data = XdfList[[k]])
  gam15 <- gam(yhat ~ s(x1, x5) + s(x2) + s(x3) + s(x4) + s(x6), data = XdfList[[k]])
  gam16 <- gam(yhat ~ s(x1, x6) + s(x2) + s(x3) + s(x4) + s(x5), data = XdfList[[k]])
  
  gam23 <- gam(yhat ~ s(x2, x3) + s(x1) + s(x4) + s(x5) + s(x6), data = XdfList[[k]])
  gam24 <- gam(yhat ~ s(x2, x4) + s(x1) + s(x3) + s(x5) + s(x6), data = XdfList[[k]])
  gam25 <- gam(yhat ~ s(x2, x5) + s(x1) + s(x3) + s(x4) + s(x6), data = XdfList[[k]])
  gam26 <- gam(yhat ~ s(x2, x6) + s(x1) + s(x3) + s(x4) + s(x5), data = XdfList[[k]])
  
  gam34 <- gam(yhat ~ s(x3, x4) + s(x1) + s(x2) + s(x5) + s(x6), data = XdfList[[k]])
  gam35 <- gam(yhat ~ s(x3, x5) + s(x1) + s(x2) + s(x4) + s(x6), data = XdfList[[k]])
  gam36 <- gam(yhat ~ s(x3, x6) + s(x1) + s(x2) + s(x4) + s(x5), data = XdfList[[k]])
  
  gam45 <- gam(yhat ~ s(x4, x5) + s(x1) + s(x2) + s(x3) + s(x6), data = XdfList[[k]])
  gam46 <- gam(yhat ~ s(x4, x6) + s(x1) + s(x2) + s(x3) + s(x5), data = XdfList[[k]])
  
  gam56 <- gam(yhat ~ s(x5, x6) + s(x1) + s(x2) + s(x3) + s(x4), data = XdfList[[k]])


  ## gamList[[k]] <- gamListK

  rsqDfK <- ## rsqDfList[[k]] <-
    rbind(
              data.frame(sim = k, rsq = rsqGamma(gamListK[[1]]$fittedValues, yhatmat), summary = "Additive", itx = "N/A"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[2]]$fittedValues, yhatmat), summary = "Additive with (x1, x2) itx", itx = "(x1, x2)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[3]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x1, x3)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[4]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x1, x4)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[5]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x1, x5)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[6]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x1, x6)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[7]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x2, x3)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[8]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x2, x4)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[9]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x2, x5)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[10]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x2, x6)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[11]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x3, x4)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[12]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x3, x5)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[13]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x3, x6)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[14]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x4, x5)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[15]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x4, x6)"),
    data.frame(sim = k, rsq = rsqGamma(gamListK[[16]]$fittedValues, yhatmat), summary = "Additive with itx (not x1 & x2)", itx = "(x5, x6)")) 

  rsqDfSummaryList[[k]] <- rsqDfK %>%
    group_by(sim, summary, itx) %>%
    summarize(rsq_mid = mean(rsq),
              rsq_lo = quantile(rsq, 0.05),
              rsq_hi = quantile(rsq, 0.95))

  

  cat(sprintf("Simulation %i out of %i...\n", k, nsims))

}


rsqDfSummary <- rsqDfSummaryList %>% plyr::rbind.fill(rsqDfSummaryList)


save.image("R-output/itx-search-big.Rdata")

load("R-output/itx-search-big.Rdata")

## Rank simulations by R^2_\gamma of (x1, x2) interaction for easier plotting later
rsqRank <- rsqDfSummary %>%
  filter(summary == "Additive with (x1, x2) itx") %>%
  dplyr::select(rsq_mid, sim) %>%
  arrange(desc(rsq_mid)) %>% 
  mutate(sim2 = 1:n()) %>%
  dplyr::select(-rsq_mid)


rsqDfSummary <- rsqDfSummary %>%
  left_join(rsqRank)

## Plot of R^2_\gamma
rsqSummaryPlot <- rsqDfSummary %>%
  ggplot() +
  geom_point(aes(sim2, rsq_mid, col = summary), 
             alpha = 0.5, position = "jitter", width = 0.2) +
  ggthemes::scale_color_colorblind() +
  labs(y = TeX("$R^2_{\\gamma}$ "), x = TeX("Simulation index \\[arranged by $R^2_{\\gamma}$ of (x1, x2)\\]")) +
  theme(legend.position = "top") +
  coord_flip()

rsqSummaryPlot

ggsave("figures/revision/rsq-point-big.png", rsqSummaryPlot,
       width = 8, height = 5, units = "in")

ggsave("figures/revision/rsq-point-big.pdf", rsqSummaryPlot, device = cairo_pdf,
       width = 8, height = 6, units = "in")


## Count the number of times that the x1x2 interaction is most
## significant, as measured by increase in R^2

for (i in 1:nsims) {

  if (i == 1) {
    wincountMid <- 0
    wincountLo <- 0
  }

  rsqx1x2Mid <- rsqDfSummary %>%
    filter(sim == i) %>%
    filter(summary == "Additive with (x1, x2) itx") %>%
    pull(rsq_mid)

  rsqx1x2Lo <- rsqDfSummary %>%
    filter(sim == i) %>%
    filter(summary == "Additive with (x1, x2) itx") %>%
    pull(rsq_lo)

  rsqOtherMid <- rsqDfSummary %>%
    filter(sim == i) %>%
    filter(summary != "Additive with (x1, x2) itx") %>%
    pull(rsq_mid) %>%
    max()

  rsqOtherHi <- rsqDfSummary %>%
    filter(sim == i) %>%
    filter(summary != "Additive with (x1, x2) itx") %>%
    pull(rsq_hi) %>%
    max()

  wincountMid <- wincountMid + (rsqx1x2Mid > rsqOtherMid)
  wincountLo <- wincountLo + (rsqx1x2Lo > rsqOtherHi)

  if (i %% 100 == 0) {
    cat(i)
    cat("\n")
  }

}

winprobMid <- wincountMid / nsims

winprobMid
sqrt(winprobMid * (1 - winprobMid) / nsims)


winprobLo <- wincountLo / nsims

winprobLo
sqrt(winprobLo * (1 - winprobLo) / nsims)

