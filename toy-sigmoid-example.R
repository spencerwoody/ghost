
library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)
library(mgcv)
library(pracma)
library(stringr)
library(mvtnorm)
library(Matrix)
library(GPFDA)
library(gridExtra)
theme_set(theme_bw(base_size = 14))

RColorBrewer::display.brewer.all()

source(here::here("R/rsq.R"))
source(here::here("R/phiGamma.R"))
source(here::here("R/rsqGamma.R"))

alpha <- 0.05

###############################################################################
                                        #           Define functions          #
###############################################################################

flogit <- function(x) log(x / (1 - x))
ilogit <- function(x) 1 / (1 + exp(-x))
ilogits <- function(x) ilogit(x * scale)

###############################################################################
                                        #             Create data             #
###############################################################################

## Length and limits of grid
M <- 50
xmax <- 2

## Create grids
x1 <- seq(-xmax, xmax, length.out = M)
x2 <- seq(-xmax, xmax, length.out = M)

sigmoidfun <- function(x1, x2) {
  x12 <- 1 * x1 + 1 * x2
  x22 <- 0.5 * x1 - 2 * x2
  1 * ilogits(x12) + 1 * ilogits(x22)
}

scale <- 2
sigma <- 0.5

## Create data
GPdf <- expand.grid(x1 = x1, x2 = x2) %>%
  mutate(x12 = 1 * x1 + 1 * x2,
         x22 = 0.5 * x1 - 2 * x2,
         f = sigmoidfun(x1, x2)) %>%
  mutate(f = f - mean(f),
         y = f + rnorm(n(), sd = sigma),
         y = y - mean(y)) 

X <- GPdf %>% dplyr::select(x1, x2) %>% as.matrix()
N <- nrow(X)

## Plot y values
ggplot(GPdf) +
   geom_raster(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette = "RdBu") +
  ## scale_fill_viridis_c() +
  coord_equal() +
  labs(title = TeX(sprintf("$\\sigma = $%1.3f", sigma)))

f <- GPdf %>% pull(f)
y <- GPdf %>% pull(y)

###############################################################################
                                        #      Fit a GP and get posterior     #
###############################################################################

GPfit <- gpr(X, y, Cov = "pow.ex")

## Empirical Bayes covariance matrix

hyper <- GPfit$hyper
fhat <- GPfit$fitted.mean
sigma2hat <- exp(GPfit$hyper$vv)

K <- GPFDA::cov.pow.ex(hyper = (hyper), Data = GPfit$train.x)
Ke <- eigen(K)

KeVec <- Ke$vectors
KeVal <- Ke$values

## How to get estimates...
fhatest <- K %*% tcrossprod(multiplyColumns(KeVec, evals), KeVec) %*% y
mean((fhat - fhatest)^2)
cor(fhat, fhatest)

## Negative eigenvalues...
fhatcovhalf <- KeVec %*% Diagonal(N, 1 / sqrt(1 / KeVal + 1 / sigma2hat))
fhatcov <- tcrossprod(fhatcovhalf)
fhatcov <- KeVec %*% Diagonal(N, 1 / (1 / KeVal + 1 / sigma2hat)) %*% t(KeVec)

NMC <- 1000
nburn <- 100
nruns <- NMC + nburn + 1

## Gibbs sampler
sigma2Samples <- rep(NA, nruns)
fhatSamples <- matrix(nrow = N, ncol = nruns)

sigma2Samples[1] <- sigma2hat
fhatSamples[, 1] <- fhat

prog <- progress_estimated(nruns)
for (i in 2:nruns) {

  ## Sample fhat
  evals <- 1 / (KeVal + sigma2Samples[i - 1])
  fhatmat <- K %*% tcrossprod(multiplyColumns(KeVec, evals), KeVec)

  fhatmean <- fhatmat %*% y
  fhatcov <- sigma2Samples[i - 1] * fhatmat

  fhatSamples[, i] <- rmvnorm(1, fhatmean, fhatcov)

  ## Sample sigma2
  RSSi <- sum((y - fhatSamples[, i - 1])^2)

  sigma2Samples[i] <- 1 / rgamma(1, N / 2, RSSi / 2)

  prog$print()$tick()
}

## Remove burn in
removeIdx <- 1:(nburn + 1)

sigma2Samples <- sigma2Samples[-removeIdx]
fhatSamples <- fhatSamples[, -removeIdx]

## save(list = ls(), file = "R-output/sigmoid-gpfit-posterior.Rdata")

## load("R-output/sigmoid-gpfit-posterior.Rdata")

## 95% credible band
ylo <- apply(fhatSamples, 1, quantile, probs = 0.025)
yhi <- apply(fhatSamples, 1, quantile, probs = 0.975)

## fhat <- GAMfit1$fitted.values
fhat <- GPfit$fitted.mean
GPdfc <- GPdf %>% mutate(fhat = fhat)

GPdf3 <- GPdfc %>%
  tidyr::pivot_longer(cols = c(f, fhat))

GPdf %>%
  mutate(fhat = fhat) %>%
  tidyr::pivot_longer(cols = c(f, fhat)) %>% 
  ggplot() +
  geom_raster(aes(x1, x2, fill = value)) +
  ## geom_path(aes(x = , y = a)) + 
  ## scale_fill_gradient2() +
  scale_fill_gradient2(low = "darkred", high = "darkblue") +
  ## scale_fill_distiller("", palette = "Spectral") +
  facet_wrap(~name) + 
  ## scale_fill_viridis_c() +
  coord_equal() +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 15))


###############################################################################
                                        #           Create summaries          #
###############################################################################

GPdf <- GPdf %>% mutate(fhat = fhat)

## Summaries, trained on fhat and y
GAMfit <- gam(fhat ~ s(x1) + s(x2), data = GPdf)
LINfit <- gam(fhat ~ (x1) + (x2), data = GPdf)

## Posterior for GAM summary
Xs <- model.matrix(GAMfit)
V <- vcov(GAMfit, dispersion = 1) 
Q <- crossprod(V, crossprod(Xs, fhatSamples))
gammaGAM <- Xs %*% Q

gammaGAM1 <- Xs[, 2:10] %*% Q[2:10, ]
gammaGAM2 <- Xs[, 11:19] %*% Q[11:19, ]

GAM1mid <- rowMeans(gammaGAM1)
GAM1lo <- apply(gammaGAM1, 1, quantile, probs = alpha / 2)
GAM1hi <- apply(gammaGAM1, 1, quantile, probs = 1 - alpha / 2)

GAM2mid <- rowMeans(gammaGAM2)
GAM2lo <- apply(gammaGAM2, 1, quantile, probs = alpha / 2)
GAM2hi <- apply(gammaGAM2, 1, quantile, probs = 1 - alpha / 2)

## Posterior for linear summary
Xsl <- model.matrix(LINfit)
Vl <- vcov(LINfit, dispersion = 1) 
Ql <- crossprod(Vl, crossprod(Xsl, fhatSamples))
gammaLIN <- Xsl %*% Ql

gammaLIN1 <- matrix(Xsl[, 2], ncol = 1) %*% Ql[2, ]
gammaLIN2 <- matrix(Xsl[, 3], ncol = 1) %*% Ql[3, ]

LIN1mid <- rowMeans(gammaLIN1)
LIN1lo <- apply(gammaLIN1, 1, quantile, probs = alpha / 2)
LIN1hi <- apply(gammaLIN1, 1, quantile, probs = 1 - alpha / 2)

LIN2mid <- rowMeans(gammaLIN2)
LIN2lo <- apply(gammaLIN2, 1, quantile, probs = alpha / 2)
LIN2hi <- apply(gammaLIN2, 1, quantile, probs = 1 - alpha / 2)

## Concatenate posterior description of summary
summaryDf <- rbind(
  data.frame(x = X[, 1], f = LIN1mid, part = "mid", term = "x[1]", model = "Linear"),
  data.frame(x = X[, 1], f = LIN1lo, part = "lo", term = "x[1]", model = "Linear"),
  data.frame(x = X[, 1], f = LIN1hi, part = "hi", term = "x[1]", model = "Linear"),
  data.frame(x = X[, 2], f = LIN2mid, part = "mid", term = "x[2]", model = "Linear"),
  data.frame(x = X[, 2], f = LIN2lo, part = "lo", term = "x[2]", model = "Linear"),
  data.frame(x = X[, 2], f = LIN2hi, part = "hi", term = "x[2]", model = "Linear"),
  data.frame(x = X[, 1], f = GAM1mid, part = "mid", term = "x[1]", model = "Additive"),
  data.frame(x = X[, 1], f = GAM1lo, part = "lo", term = "x[1]", model = "Additive"),
  data.frame(x = X[, 1], f = GAM1hi, part = "hi", term = "x[1]", model = "Additive"),
  data.frame(x = X[, 2], f = GAM2mid, part = "mid", term = "x[2]", model = "Additive"),
  data.frame(x = X[, 2], f = GAM2lo, part = "lo", term = "x[2]", model = "Additive"),
  data.frame(x = X[, 2], f = GAM2hi, part = "hi", term = "x[2]", model = "Additive")
)

## Plot it
summaryPlot <- summaryDf %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_line(aes(x, f, col = model, lty = part), size=0.75, alpha = 0.9) +
  facet_wrap(~term, labeller=label_parsed) +
  scale_linetype_manual(values = c(1, 2, 2), guide = FALSE) +
  scale_color_brewer("Summary", palette = "Set2") + 
  labs(title = "Partial effects from summaries",
       x = TeX("x_j"), y = TeX("$f_j(x_j)$"),
       tag = "(b)") +
  theme(legend.position = "bottom")


summaryPlot

ggsave("figures/sigmoid/01sigmoid-univariate-plots.pdf", summaryPlot,
       width = 4, height = 10.5)


summaryPlot2 <- summaryDf %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_line(aes(x, f, col = model, lty = part), size=0.75, alpha = 0.9) +
  facet_wrap(~term, labeller=label_parsed, ncol = 1) +
  scale_linetype_manual(values = c(1, 2, 2), guide = FALSE) +
  scale_color_brewer("Summary", palette = "Set2") + 
  labs(title = "Summary partial effects",
       x = TeX("x_j"), y = TeX("$f_j(x_j)$"),
       tag = "(b)") +
  theme(legend.position = "bottom")

summaryPlot2

###############################################################################
                                        #        Compare bivariate fits       #
###############################################################################

GPdf <- GPdf %>%
  mutate(fhat = fhat,
         fhatgam = GAMfit$fitted.values,
         fhatlin = LINfit$fitted.values)

GPdf2 <- GPdf %>%
  tidyr::pivot_longer(cols = c(f, y, fhat, fhatgam, fhatlin))

GPdf2 <- GPdf2 %>%
  mutate(name = case_when(
           name == "f" ~ "True function",
           name == "fhat" ~ "Estimated function",
           name == "y" ~ "Observations",
           name == "fhatgam" ~ "Additive summary",
           name == "fhatlin" ~ "Linear summary",
           TRUE ~ name
         ),
         name = factor(name, levels = c(
                               "True function", "Observations",
                               "Estimated function", "Additive summary",
                               "Linear summary"
                             )))

compPlot1 <- GPdf2 %>% filter(name!="f2") %>% 
  ggplot() +
  geom_raster(aes(x1, x2, fill = value)) +
  ## scale_fill_gradient2(low = "darkred", high = "darkblue") +
  scale_fill_gradient2(TeX("$f(x_1, x_2)$"), low = "#660000", high = "#1D2951") +
  ## scale_fill_distiller(palette = "Spectral") + 
  facet_wrap(~name, ncol = 1) + 
  ## scale_fill_viridis_c() +
  coord_equal() +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 15)) +
  labs(title = "Estimated function\nand summaries",
       x = TeX("$x_1$"), y = TeX("$x_2$"),
       subtitle = TeX(sprintf("$\\sigma^2 = $%1.3f", sigma^2)),
       tag = "(a)") +
  theme(legend.position = "left")

compPlot1



ggsave("figures/sigmoid/sigmoid01-estimated-functions.pdf", compPlot1,
       width = 4, height = 10.5)



compPlot12 <- GPdf2 %>% filter(name!="f2") %>% 
  ggplot() +
  geom_raster(aes(x1, x2, fill = value)) +
  ## scale_fill_gradient2(low = "darkred", high = "darkblue") +
  scale_fill_gradient2(TeX("$f(x_1, x_2)$"), low = "#660000", high = "#1D2951") +
  ## scale_fill_distiller(palette = "Spectral") + 
  facet_wrap(~name, nrow = 1) + 
  ## scale_fill_viridis_c() +
  coord_equal() +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 8)) +
  labs(title = "Estimated function and summaries",
       x = TeX("$x_1$"), y = TeX("$x_2$"),
       subtitle = TeX(sprintf("$\\sigma^2 = $%1.3f", sigma^2)),
       tag = "(a)") +
  theme(legend.position = "right")

compPlot12


###############################################################################
                                        #         Summary diagnostics         #
###############################################################################

## Create summary statistics

rsqGAM <- rsqGamma(gammaGAM, fhatSamples)
rsqLIN <- rsqGamma(gammaLIN, fhatSamples)

phiGAM <- phiGamma(y, gammaGAM, sqrt(sigma2Samples))
phiLIN <- phiGamma(y, gammaLIN, sqrt(sigma2Samples))

mean(rsqGAM)
mean(rsqLIN)

mean(phiGAM)
mean(phiLIN)

ggplot() +
  geom_density(aes(rsqGAM, col = "GAM")) +
  geom_density(aes(rsqLIN, col = "LIN"))

ggplot() +
  geom_density(aes(phiGAM, col = "GAM")) +
  geom_density(aes(phiLIN, col = "LIN"))

###############################################################################
                                        #         Partial derivatives         #
###############################################################################

Q <- GPfit$Q
Qinv <- solve(Q)

mean((Q - K - diag(sigma2hat, N))^2)

predGPfit <- function(x) {
  ## gppredict(GPfit, Data.new = data.frame(x1 = x[1], x2 = x[2]))
  Kjoint <- cov.pow.ex(hyper,
                       Data = GPfit$train.x,
                       Data.new = mypoint <- matrix(c(x1 = x[1], x2 = x[2]), nrow = 1))
  as.numeric(Kjoint %*% Qinv %*% y)
}

predGAMfit <- function(x) {
  predict(GAMfit, newdata=data.frame(x1 = x[1], x2 = x[2]))
}

predsigmoidfun <- function(x) {
  sigmoidfun(x[1], x[2])
}

GPdf <- GPdf %>%
  mutate(gradx1f = NA,
         gradx1fhat = NA,
         gradx1gam = NA,
         gradx1lin = NA,
         gradx2f = NA,
         gradx2fhat = NA,
         gradx2gam = NA,
         gradx2lin = NA)

prog <- progress_estimated(nrow(GPdf))
for (i in 1:nrow(GPdf)) {

  x1i <- GPdf$x1[i]; x2i <- GPdf$x2[i]

  gradsigmoidfun <- grad(predsigmoidfun, c(x1i, x2i))
  
  GPdf$gradx1f[i] = gradsigmoidfun[1]
  GPdf$gradx2f[i] = gradsigmoidfun[2]

  gradGP <- grad(predGPfit, c(x1i, x2i))

  GPdf$gradx1fhat[i] <- gradGP[1]
  GPdf$gradx2fhat[i] <- gradGP[2]

  gradGAM <- grad(predGAMfit, c(x1i, x2i))

  GPdf$gradx1gam[i] <- gradGAM[1]
  GPdf$gradx2gam[i] <- gradGAM[2]

  GPdf$gradx1lin[i] = coef(LINfit)[2]
  GPdf$gradx2lin[i] = coef(LINfit)[3]
 
  prog$tick()$print()
}

## save(list=ls(), file = "R-output/sigmoid-example-posterior-grad.Rdata")
load("R-output/sigmoid-example-posterior-grad.Rdata")

## Make this a "long" dataframe, for plotting purposes
gradDf <- GPdf %>%
  tidyr::pivot_longer(cols = c(gradx1f, gradx2f,
                               gradx1fhat, gradx2fhat,
                               gradx1gam, gradx2gam,
                               gradx1lin, gradx2lin)) %>%
  mutate(term = case_when(
           str_detect(name, "x1") ~ "x1",
           TRUE ~ "x2"
         ),
        term2 = case_when(
           str_detect(name, "x1") ~ "1",
           TRUE ~ "2"),
         model = case_when(
           str_detect(name, "fhat") ~ "Estimated function",
           str_detect(name, "gam") ~ "Additive summary",
           str_detect(name, "lin") ~ "Linear summary",
           TRUE ~ "True function"
         ),
         model = factor(model, levels = c(
                                 "True function", "Estimated function",
                                 "Additive summary", "Linear summary"
                               )))

## Gradient plot
gradPlot <- gradDf %>%
  ggplot() +
  geom_raster(aes(x1, x2, fill = value)) +
  facet_grid(cols = vars(term2), rows = vars(model),
             labeller=label_bquote(cols = partialdiff*f/partialdiff*x[.(term2)])) +  
  scale_fill_gradient2(name = TeX("$\\partial f / \\partial x_j$"),
                       low = "darkorchid4", high = "darkgreen") +
  coord_equal() +
  guides(fill = guide_colourbar(barheight = 1, barwidth = 15)) +
  labs(title = "True and estimated\npartial derivatives from summaries",
       x = TeX("$x_1$"), y = TeX("$x_2$"),
       tag = "(c)") +
  theme(legend.position = "top")

gradPlot

ggsave("figures/sigmoid/sigmoid01-partial-derivatives.pdf", gradPlot,
       width = 8, height = 9)

## Move the legend
gradPlot2 <- gradDf %>%
  ggplot() +
  geom_raster(aes(x1, x2, fill = value)) +
  facet_grid(rows = vars(term2), cols = vars(model),
             labeller=label_bquote(rows = partialdiff*f/partialdiff*x[.(term2)])) +  
  scale_fill_gradient2(name = TeX("$\\partial f / \\partial x_j$"),
                       low = "darkorchid4", high = "darkgreen") +
  coord_equal() +
  guides(fill = guide_colourbar(barheight = 15)) +
  labs(title = "True and estimated partial derivatives from summaries",
       x = TeX("$x_1$"), y = TeX("$x_2$"),
       tag = "(c)") +
  theme(legend.position = "right")

gradPlot2

 ###############################################################################
                                        #        Concatenate the plots        #
###############################################################################

lay <- rbind(
  c(1, 1, 2, 2),
  c(1, 1, 3, 3),
  c(1, 1, 3, 3)
)

grobList <- list(compPlot1, summaryPlot, gradPlot)

grobCat <- arrangeGrob(compPlot1, summaryPlot, gradPlot, layout_matrix = lay)

grid.arrange(grobCat)

ggsave("figures/sigmoid01-concatenated.pdf", grobCat,
       width = 10.3, height = 11.9, units =  "in")

## landscape mode

lay2 <- rbind(
  c(1, 1, 1, 1, 1, 1, 1),
  c(1, 1, 1, 1, 1, 1, 1),
  c(1, 1, 1, 1, 1, 1, 1),
  c(2, 2, 3, 3, 3, 3, 3),
  c(2, 2, 3, 3, 3, 3, 3),
  
  c(2, 2, 3, 3, 3, 3, 3),
  c(2, 2, 3, 3, 3, 3, 3)
)

grobCat2 <- arrangeGrob(compPlot12, summaryPlot2, gradPlot2, layout_matrix = lay2)

grid.arrange(grobCat2)

## ggsave("figures/sigmoid-example-stacked2.pdf", grobCat2)

ggsave("figures/sigmoid-example-stacked2.pdf", grobCat2,
       width = 13, height=9, units = "in")

## END
