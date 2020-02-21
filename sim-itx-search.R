
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
library(gridExtra)
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
                                        #             simulations for GP      #
###############################################################################

nsims <- 21

XList <- vector("list", nsims)
fList <- vector("list", nsims)
yList <- vector("list", nsims)
GPfitList <- vector("list", nsims)
postList <- vector("list", nsims)

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

  rsqDfK <- rsqDfList[[k]] <- rbind(
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

## Concatenate lists of dataframes into dataframes
rsqDf <- rsqDfList %>% plyr::rbind.fill(rsqDfList)

rsqDfSummary <- rsqDfSummaryList %>% plyr::rbind.fill(rsqDfSummaryList)

## Rank the simulations by R^2 of "correct" interaction (but put
## "last" simulation first, since this one is used for plotting)
rsqRank <- rsqDfSummary %>%
  filter(summary == "Additive with (x1, x2) itx") %>%
  dplyr::select(rsq_mid, sim) %>%
  mutate(isLast = (sim == nsims)  * 1) %>%
  arrange(desc(isLast), desc(rsq_mid)) %>% 
  mutate(sim2 = 1:n()) %>%
  dplyr::select(-rsq_mid)

rsqDf <- rsqDf %>% left_join(rsqRank) 

rsqDfSummary <- rsqDfSummary %>% left_join(rsqRank) 

save.image("R-output/itx-search.Rdata")

load("R-output/itx-search.Rdata")


###############################################################################
                                        #  Plot results of sims of summaries  #
###############################################################################



## R^2_gamma comparison for one simulation
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

rsqPlot1Facet <- rsqDfK %>%
  mutate(summary = factor(summary,
                          levels = c("Additive",
                                     "Additive with itx (not x1 & x2)",
                                     "Additive with (x1, x2) itx"))) %>% 
  filter(!(str_detect(itx, "x4, ") | str_detect(itx, "x5, "))) %>%
  ggplot() +
  geom_density(aes(rsq, col = summary, group = itx), size = 0.75) +
  facet_wrap(~summary, ncol = 1, scales = "free_y") + 
  labs(x = TeX("$R^2_\\gamma$")) +
  ## ggthemes::scale_color_colorblind() +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[3], cbbPalette[2])) +
  theme(legend.position = "none") 

rsqPlot1Facet

ggsave("figures/revision/rsqplot1-facet.pdf", rsqPlot1Facet,
       width = 6.5, height = 5.25, units = "in")

## R^2_gamma comparison for all simulations
rsqRidge <- rsqDf %>%
  filter(!is.na(sim2)) %>%
  ## left_join(rsqRank) %>%
  ggplot() +
  geom_density_ridges(aes(x = rsq, y = as.factor(sim2),
                          col = itx, fill = summary), alpha = 0.5) +
  scale_color_manual(values = rep("black", 50)) +
  ggthemes::scale_fill_colorblind() +
  theme_cowplot(font_size = 18) + 
  theme(legend.position = "top") +
  guides(color = "none") +
  labs(x = TeX("$\\R^2_{\\gamma}$"),
       y = TeX("simulation index \\[arranged by $R^2_{\\gamma}$ of (x1, x2)\\]"))

rsqRidge

ggsave("figures/revision/rsq-ridges.pdf", rsqRidge, width = 10, height = 8, units = "in")



#######################################################################################
                                        # Plot the additive and
                                        # interaction summaries
#######################################################################################

## Purely additive plot
gamListK[[1]]$gamPlot

ggsave("figures/revision/gam1-plot.pdf", gamListK[[1]]$gamPlot +ylab(TeX("$h_j(x_j)$")),
       width = 7, height = 5, units = "in")

## Summary after adding (x1, x2) interaction: only the additive part
gamplot21 <- gamListK[[2]]$gamDf %>%
  filter(!(term %in% c("x1", "x2"))) %>%
  ggplot() +
  geom_hline(yintercept = 0) + 
  geom_ribbon(aes(x_j, ymin = fx_j_lo, ymax = fx_j_hi),
              fill = "grey50", alpha = 0.5) +
  geom_line(aes(x_j, fx_j_mean), col = "firebrick3") +
  geom_rug(aes(x_j, fx_j_mean), sides = "b", alpha = 0.25) + 
  facet_wrap(~term, scale = "free_x") +
                                        # theme_minimal(base_size = 14) + 
  theme_half_open(font_size = 16) + 
  labs(x = TeX("$x_j$"), y = TeX("$f_j(x_j)$")) 

gamplot21

ggsave("figures/revision/gam2-plot-additive.pdf", gamplot21,
       width = 6, height = 4.5, units = "in")

## Plot the interaction surface

## Resolution of image
seqlen <- 64

## Create a grid of X1, X2
X1seq <- seq(min(XList[[k]][, 1]), max(XList[[k]][, 1]), length.out = seqlen)
X2seq <- seq(min(XList[[k]][, 2]), max(XList[[k]][, 2]), length.out = seqlen)

X1X2 <- expand.grid(X1seq, X2seq)

Xnew <- data.frame(
  x1 = X1X2[, 1],
  x2 = X1X2[, 2],
  x3 = 0, x4 = 0, x5 = 0, x6 = 0
)


## Model matrix from fit
Xgam12 <- model.matrix(gam12)

## Model matrix for new data
Xgam12new <- model.matrix(gam12, newdata = Xnew)

## Posterior for coefficients
Vgam12 <- vcov(gam12, dispersion = 1)
Qgam12 <- crossprod(Vgam12, crossprod(Xgam12, yhatmat))

## Posterior for fit
itxtermfit <- Xgam12new[, 2:30] %*% Qgam12[2:30, ]

## Data frame for plotting
itxplotdf <- data.frame(
  x1 = X1X2[, 1],
  x2 = X1X2[, 2],
  Mean = rowMeans(itxtermfit),
  Var = apply(itxtermfit, 1, var)
)

## Plot posterior mean and variance
fx1x2Plot <- ggplot() +
  geom_point(data = NULL, aes(XList[[k]][, 1], XList[[k]][, 2]), alpha = 0.5) + 
  geom_raster(data = itxplotdf, aes(x1, x2, fill = Mean), alpha = 0.8) +
  scale_fill_viridis_c(TeX("$h_{1,2}(x_1, x_2)$")) +
  ## scale_fill_gradient2(TeX("$f(x_1, x_2)$")) +
  theme_half_open(font_size=16) + 
  theme(legend.position = "top") +
  guides(fill = guide_colourbar(barwidth = 12, barheight = 0.75, title.position = "top")) +
  coord_equal() +
  labs(x = "x1", y = "x2")

fx1x2Plot


fx1x2VarPlot <- ggplot() +
  geom_point(data = NULL, aes(XList[[k]][, 1], XList[[k]][, 2]), alpha = 0.5) + 
  geom_raster(data = itxplotdf, aes(x1, x2, fill = Var), alpha = 0.8) +
  scale_fill_viridis_c(TeX("$var(h_{1,2}(x_1, x_2) | y)$"), option = "C") +
  theme_half_open(font_size=16) + 
  theme(legend.position = "top") +
  guides(fill = guide_colourbar(barwidth = 12, barheight = 0.75, title.position = "top")) +
  coord_equal() +
  labs(x = "x1", y = "x2")

fx1x2VarPlot

## Arrange two plots of partially additive
mygrob <- arrangeGrob(fx1x2Plot, fx1x2VarPlot, nrow = 1)

grid.arrange(mygrob)

ggsave("figures/revision/gam2-plot-bivariate.pdf", mygrob, device = cairo_pdf,
       width = 8, height = 4.5)


## Save the regression tree
rpart.plot::rpart.plot(treeList[[k]])

pdf("figures/revision/gam-tree-1.pdf", width = 4.5, height = 3)
rpart.plot::rpart.plot(treeList[[k]], tweak = 1.2)
dev.off()

