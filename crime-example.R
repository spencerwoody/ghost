
library(here)
library(horseshoe)
library(dplyr)
library(MASS)
library(ggplot2)
library(lars)
library(latex2exp)
library(reshape2)
library(mvtnorm)
theme_set(theme_bw(base_size = 16))

source(here("R/sparseLinearProj.R"))
source(here("R/sparseLinearFuns.R"))
source(here("R/rho.R"))
source(here("R/phiGamma.R"))
source(here("R/rsqGamma.R"))
source(here("R/phiGammaLinear.R"))
source(here("R/rsqGammaLinear.R"))
source(here("R/psi.R"))

## Confidence level for credible intervals
alpha <- 0.05

###############################################################################
                                        #              Data prep              #
###############################################################################

data(UScrime)

## Response variable, scaled
y <- UScrime %>%
  pull(y) %>%
  scale()

## Covariates
X <- UScrime[, -which(colnames(UScrime) == "y")] %>%
  as.matrix()

## log-transform and scale the data
X[, -2] <- log(X[, -2])

X <- scale(X)

## Variable names
varnames <- colnames(X)

varnamesDf <- data.frame(
  Var1 = seq_along(varnames),
  varname = varnames
)

varnamesFactor <- varnames %>% factor(levels = varnames)

###############################################################################
                                        #       Fit model with horseshoe      #
###############################################################################

## Fit Bayesian model
fit <- horseshoe(y, X, method.tau="halfCauchy",
                 method.sigma="Jeffreys", nmc=45000, thin=5)

## Posterior samples of beta; each covariate is one row
betaSamples <- fit$BetaSamples

## Posterior mean
betaBar <- rowMeans(betaSamples)

## Estimated variance of residuals
fit$Sigma2Hat

sigma2Samples <- fit$Sigma2Samples

## Predicted value of y
yFit <- rowMeans(X %*% betaSamples)

## Dataframe of beta samples
betaSamplesMelt <- betaSamples %>%
      melt() %>%
      dplyr::select(-Var2) %>%
      left_join(varnamesDf, by = "Var1") %>%
      dplyr::select(-Var1)

betaSamplesSummary <- betaSamplesMelt %>%
  group_by(varname) %>%
  summarize(mid = mean(value),
            lo = quantile(value, alpha / 2),
            hi = quantile(value, 1 - alpha / 2))


###############################################################################
                                        #   Perform summary search and proj.  #
###############################################################################

sparesProj <- sparseLinearProj(X, y,
                               betaSamples, sigma2Samples,
                               varnames = varnames)

sparesProj$betaProjSummary %>% glimpse()

## Change point estimate to DSS estimate
for (j in 1:ncol(sparesProj$betaLambda)) {
  for (k in 1:length(varnames)) {
    sparesProj$betaProjSummary <- sparesProj$betaProjSummary %>%
      mutate(mid = case_when(
               as.character(varname) == varnames[k] &
               as.character(modelSize) == as.character(j) ~ sparesProj$betaLambda[k, j],
               TRUE ~ mid
             ))
  }
}

names(sparesProj)

###############################################################################
                                        #       All projected posteriors      #
###############################################################################

allProjPostPlot <- ggplot() +
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_violin(data = sparesProj$betaProjDf %>%
                mutate(modelSize = factor(modelSize,
                                          levels = rev(levels(modelSize)))),
              aes(modelSize, value), col = "grey50", scale = "width") +
  geom_errorbar(data = sparesProj$betaProjSummary %>% filter(mid != 0),
                aes(modelSize, ymin = lo, ymax = hi), width = 0.4) +
  geom_point(data = sparesProj$betaProjSummary %>% filter(mid != 0),
             aes(modelSize, mid, col = "Mean of proj. posterior")) +
  geom_point(data = sparesProj$betaLambdaDf %>% filter(value != 0),
             aes(modelSize, value, col = "Point estimate")) +
  facet_wrap(~varname) +
  labs(title = "Projected posteriors",
       x = "No. of variables in summary", y = TeX("$\\beta_j$")) +
  scale_color_brewer("", palette = "Set1") + 
  theme(legend.position = "bottom") + 
  coord_flip()

allProjPostPlot

ggsave("figures/crime-projected-posteriors-flip.pdf", allProjPostPlot,
       width = 13.5, height = 13, units = "in")

###############################################################################
                                        #               Po1, Po2              #
###############################################################################

projDf <- sparesProj$betaProjDf %>%
  mutate(modelSize = factor(modelSize,
                            levels = rev(levels(modelSize)))) %>%
  filter(!(modelSize %in% c("1", "2", "3", "4", "5"))) %>%
  filter(varname %in% c("Po1", "Po2")) %>%
  filter(value!=0)

projSummary <- sparesProj$betaProjSummary %>% filter(mid != 0) %>%
  filter(!(modelSize %in% c("1", "2", "3", "4", "5"))) %>%
  filter(varname %in% c("Po1", "Po2")) %>%
  filter(mid != 0)

PoPlot <- ggplot() +
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_violin(data = projDf, aes(modelSize, value),
              alpha = 0.75, col = "grey50", scale = "width") +
  geom_errorbar(data = projSummary, aes(modelSize, ymin = lo, ymax = hi),
                width = 0.4) + 
  geom_point(data = projSummary, aes(modelSize, mid),
             col = "red") +
  facet_wrap(~varname) +
  labs(title = "Projected posteriors for two highly collinear variables",
       x = "Model Size", y = TeX("$\\beta_j$")) +
  coord_flip()

PoPlot

ggsave("figures/crime-Po1-Po2.pdf", PoPlot, width = 12, height = 6.8, units = "in")

###############################################################################
                                        #        Summarization metrics        #
###############################################################################

## Rsq_gamma
rsqPlot <- ggplot() +
  geom_hline(yintercept = 0.95, lty = "dashed") + 
  geom_violin(data = sparesProj$summaryDfLong %>%
                filter(stat %in% c("rsq_gamma")),
              aes(modelSize, value), alpha = 0.1, col = "grey75",
              scale = "width") +
  geom_errorbar(data = sparesProj$summaryDfCI %>%
                  filter(stat %in% c("rsq_gamma")),
                aes(modelSize, ymin = lo, ymax = hi), width = 0.4) + 
  geom_point(data = sparesProj$summaryDfCI %>%
               filter(stat %in% c("rsq_gamma")),
             aes(modelSize, mid), col = "firebrick3") +
  labs(title = "Variability in full linear model explained by summary",
       x = "No. of variables in summary", y = TeX("$R^2_{\\gamma}$")) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

rsqPlot

## 11.5 x 5.25
ggsave("figures/crime-rsqgamma.pdf", rsqPlot,
       width = 11.5, height = 5.25, units = "in")

## phi_gamma
phiPlot <- ggplot() +
  geom_hline(yintercept = 0.05, lty = "dashed") + 
  geom_violin(data = sparesProj$summaryDfLong %>%
                filter(stat %in% c("phi_gamma")),
              aes(modelSize, value), alpha = 0.1, col = "grey75",
              scale = "width") +
  geom_errorbar(data = sparesProj$summaryDfCI %>%
                  filter(stat %in% c("phi_gamma")),
                aes(modelSize, ymin = lo, ymax = hi), width = 0.4) + 
  geom_point(data = sparesProj$summaryDfCI %>%
               filter(stat %in% c("phi_gamma")),
             aes(modelSize, mid), col = "firebrick3") +
  labs(title = "Inflation of residual stdev using summary",
       x = "No. of variables in summary", y = TeX("$\\phi_{\\gamma}$")) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

phiPlot

ggsave("figures/crime-phigamma.pdf", phiPlot,
       width = 11.5, height = 5.25, units = "in")

###############################################################################
                                        #       Get marginal posteriors       #
###############################################################################

selectedVars <- sparesProj$selectedVars[[6]]

Xp <- X[, sparesProj$selectedVars[[6]]]

betaMargp <- betaSamplesMelt %>%
  filter(varname %in% varnames[selectedVars])

betaMargSump <- betaSamplesSummary %>%
  filter(varname %in% varnames[selectedVars])

###############################################################################
                                        #      Refitted model, flat prior     #
###############################################################################

postMean <-  solve(crossprod(Xp), crossprod(Xp, y))

yhatp <- Xp %*% postMean

sig2hat <- var(y - yhatp)

covmat <- solve(crossprod(Xp)) * as.numeric(sig2hat)

betaSamples2 <- rmvnorm(9000, postMean, covmat) %>% t()

varnamesDf2 <- data.frame(
  Var1 = 1:length(sparesProj$selectedVars[[6]]),
  varname = varnames[sparesProj$selectedVars[[6]]]
)

betaSamples2Df <- betaSamples2 %>%
  melt() %>%
  dplyr::select(-Var2) %>%
  left_join(varnamesDf2, by = "Var1") %>%
  dplyr::select(-Var1)

mid <- apply(betaSamples2, 1, mean)
lo <- apply(betaSamples2, 1, quantile, probs = alpha / 2)
hi <- apply(betaSamples2, 1, quantile, probs = 1 - alpha / 2)

dfrefit <- data.frame(varname = varnames[sparesProj$selectedVars[[6]]],
           mid = mid, lo = lo, hi = hi) %>%
  mutate(mode = "Refitted (flat prior)")

###############################################################################
                                        #          Compare posteriors         #
###############################################################################

## Projected model
betaSamplesProj <- sparesProj$betaProjDf %>%
  filter(modelSize == sparesProj$modelSize[6],
         value!=0) %>%
  dplyr::select(-modelSize)

dfproj <- sparesProj$betaProjSummary %>% filter(modelSize == "6") %>%
  mutate(mode = "Projected") %>% ungroup() %>% 
  dplyr::select(-modelSize) %>%
  filter(mid != 0)

## Concatenate dataframes
betaSamplesCat2 <- rbind(betaSamples2Df %>% mutate(mode = "Refitted (flat prior)"),
                        betaSamplesProj %>% mutate(mode = "Projected"),
                        betaMargp %>% mutate(mode = "Marginal")) %>%
  mutate(mode = factor(mode,
                       levels = c("Marginal", "Projected",
                                  "Refitted (flat prior)")))

dfall <- rbind(dfrefit, dfproj, betaMargSump %>% mutate(mode = "Marginal")) %>%
  mutate(mode = factor(mode,
                       levels = c("Marginal", "Projected",
                                  "Refitted (flat prior)")))

## Create plot
mywidth <- 1 #Width separating violins

postCompPlot <- betaSamplesCat2 %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_violin(aes(varname, value, fill = mode), scale = "width",
              position = position_dodge(width=mywidth), alpha = 0.25) +
  geom_errorbar(data = dfall, aes(varname, ymin = lo, ymax = hi,
                                   group = mode, col = mode),
                position = position_dodge(width=mywidth),
                width = 0.3,
                size = 1) +
  geom_point(data = dfall, aes(varname, mid, col = mode),
             position = position_dodge(width=mywidth), size = 2) +
  scale_color_brewer("Posterior", palette = "Dark2") +
  scale_fill_brewer("Posterior", palette = "Dark2") + 
  labs(title = "Posterior uncertainty estimates for final selected sparse set",
       x = "Variable", y = TeX("$\\beta_j")) +
  theme(legend.position = "bottom") 

postCompPlot

ggsave("figures/crime-post-proj-refit-marg.pdf", postCompPlot,
       width = 11.6, height = 6.9, units = "in")

## END

