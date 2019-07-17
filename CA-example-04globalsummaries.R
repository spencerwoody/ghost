
library(here)
library(ggplot2)
library(dplyr)
library(rpart)
library(rpart.plot)
library(latex2exp)
library(GPFDA)
library(Matrix)
library(mvtnorm)

theme_set(theme_bw(base_size = 16))

source(here("R/rsq.R"))
source(here("R/phiGamma.R"))
source(here("R/rsqGamma.R"))

load("R-output/CA-GPfit.Rdata") 
load("R-output/GPgibbs-samples.Rdata")

rsq(y, yhat)

###############################################################################
                                        #        Create linear summary        #
###############################################################################

## Confidence level
alpha <- 0.05

## Design matrix, including intercept
xp <- cbind(1, x)

XtXinv <- solve(crossprod(cbind(xp)))

## Point estimate, and posterior of projected linear summary
betaPoint <- crossprod(XtXinv, crossprod(xp, yhat))
betaPost <- crossprod(XtXinv, crossprod(xp, yhatSamples))

## Posterior mean and credible intervals 
betaMean <- rowMeans(betaPost)
betalo <- apply(betaPost, 1, quantile, probs = alpha / 2)
betahi <- apply(betaPost, 1, quantile, probs = 1 - alpha / 2)

cor(betaMean, betaPoint)

###############################################################################
                                        #     Compare OLS and projections     #
###############################################################################

## OLS regression
ols <- lm(y~x)

coef(ols)
betaMean

## Confidence / credible intervals
ciProj <- t(apply(betaPost, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2)))
ciOLS <- confint(ols)

ROWNAMES <- rownames(ciProj)
ROWNAMES[1] <- "(Intercept)"

##  Dataframe for plotting
linDf <- rbind(
  data.frame(
    mode = "Projected",
    term = ROWNAMES,
    est = betaMean,
    lo = ciProj[, 1],
    hi = ciProj[, 2]
  ),
  data.frame(
    mode = "OLS",
    term = ROWNAMES,
    est = coef(ols),
    lo = ciOLS[, 1],
    hi = ciOLS[, 2]
  )
) %>%
  mutate(term = factor(term, levels = ROWNAMES))

linPlot <- linDf %>%
  ggplot() +
  geom_errorbar(aes(mode, ymin = lo, ymax = hi), width = 0.2) +
  geom_point(aes(mode, est), col = "firebrick3") +
  coord_flip() +
  facet_wrap(~term, scales = "free_x") +
  labs(title = "Comparing coefficients from in projected linear summary and OLS",
       subtitle = sprintf("%2.0f%% credible / confidence intervals", 100 * (1 - alpha)),
       x = "", y = "Estimate and 95% CI") 

linPlot

ggsave("figures/GP-linear-summary.pdf", linPlot,
       width = 9.9, height = 5.7, units = "in")

###############################################################################
                                        #        Fit the GAM, get posterior   #
###############################################################################

yhat <- GPfit$fitted.mean

gamFit <- gam(yhat ~
                s(Median_household_income) +
                s(POPULATION) +
                s(Median_rooms) +
                s(LONGITUDE) +
                s(LATITUDE),
              data = califDf)

## Model matrix and covariance matrix
Xs <- model.matrix(gamFit)
V <- vcov(gamFit, dispersion = 1)

## Point summary and posterior of projection
q <- crossprod(V, crossprod(Xs, yhat))
Q <- crossprod(V, crossprod(Xs, yhatSamples))

## Posterior of fitted values from projections
gammaSamples <- Xs %*% Q

gamma <- predict(gamFit, newdata = califDf, type = "response")
gammaTerms <- predict(gamFit, newdata = califDf, type = "terms",
                      se.fit = TRUE)

###############################################################################
                                        #  Make a plot looping through terms  #
###############################################################################

## Make plot
allterminds <- list(2:10, 11:19, 20:28, 29:37, 38:46)

dfList <- vector("list", length = length(allterminds)) 
gammaSElist <- vector("list", length = length(allterminds)) 

for (k in 1:length(allterminds)) {

  terminds <- allterminds[[k]]

  gammaTermSamples <- Xs[, terminds] %*% Q[terminds, ]

  gammaTermMean <- rowMeans(gammaTermSamples)
  gammaTermLo <- apply(gammaTermSamples, 1, quantile, probs = alpha / 2)
  gammaTermHi <- apply(gammaTermSamples, 1, quantile, probs = 1 - alpha / 2)

  dfList[[k]] <- data.frame(
    j = terms[k],
    xj = x[, k],
    fjxj = gammaTermMean,
    lo = gammaTermLo,
    hi = gammaTermHi
  )

  gammaSElist[[k]] <- data.frame(
    j = terms[k],
    xj = x[, k],
    fjxj = gammaTerms[[1]][, k],
    lo = gammaTerms[[1]][, k] - 1.96 * gammaTerms[[2]][, k],
    hi = gammaTerms[[1]][, k] + 1.96 * gammaTerms[[2]][, k]
  )
  
  print(k)
}

gammaSEDf <- gammaSElist %>% plyr::rbind.fill()

gammaDf <- dfList %>% plyr::rbind.fill()

## Predicted GAM and 95\% credible intervals
GAMplot <- gammaDf %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_rug(aes(xj, fjxj), sides = "b", alpha = 0.1) +
  ## geom_line(aes(xj, lo), col = "firebrick3", lty = "dashed") +
  ## geom_line(aes(xj, hi), col = "firebrick3", lty = "dashed") +
  geom_ribbon(aes(xj, ymin = lo, ymax = hi), fill = "grey60", alpha = 0.5) +
  geom_line(aes(xj, fjxj), col = "firebrick3") +
  facet_wrap(~j, scales = "free_x") +
  labs(title = "Projected additive summary of GP fit",
       subtitle = "Using posterior draws of GP",
       x = "",
       y = TeX("$f_j(x_j)$")) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) 

###############################################################################
                                        #      Plot CART of GAM residuals     #
###############################################################################

gamres <- as.numeric(yhat - gamma)

resGAMtree <- rpart(gamres ~ ., data = califDf,
                    control = rpart.control(maxdepth=5, cp = 0.0001))

resGAMtree <- rpart(gamres ~ ., data = califDf,
                    control = rpart.control(maxdepth=4, cp = 0.0001))

## resGAMtree <- rpart(gamres ~ ., data = califDf)

rpart.plot(resGAMtree,
           box.palette = "Greys",
           main = "Residuals of GAM fit to GP posterior mean",
           cex = 1)


pdf(here("figures/GP-whole-gam-res.pdf"), width = 12, height = 12)
rpart.plot(resGAMtree,
           box.palette = "Greys",
           main = "Residuals of GAM fit to GP posterior mean",
           cex = 1)
dev.off()

## ggsave(here("figures/img3/GP-gam-res-.pdf"))


###############################################################################
                                        #    Add LON,LAT interaction to GAM   #
###############################################################################

gamFit2 <- gam(yhat ~
                s(Median_household_income) +
                s(POPULATION) +
                s(Median_rooms) +
                s(LONGITUDE, LATITUDE),
              data = califDf)

gamma2 <- predict(gamFit2, newdata = califDf, type = "response")

Xs2 <- model.matrix(gamFit2)

Xs2p <- predict(gamFit2, type = "lpmatrix")

V2 <- vcov(gamFit2, dispersion = 1)
q2 <- crossprod(V2, crossprod(Xs2, yhat))
Q2 <- crossprod(V2, crossprod(Xs2, yhatSamples))

gammaSamples2 <- Xs2 %*% Q2

Xs2[1, 1:10]

dim(Xs2)

colnames(Xs2)

## Loop through terms for plotting
allterminds2 <- list(2:10, 11:19, 20:28)

dfList2 <- vector("list", length = length(allterminds2)) 
gammaSElist2 <- vector("list", length = length(allterminds2)) 

for (k in 1:length(allterminds2)) {

  terminds <- allterminds2[[k]]

  gammaTermSamples <- Xs2[, terminds] %*% Q2[terminds, ]

  gammaTermMean <- rowMeans(gammaTermSamples)
  gammaTermLo <- apply(gammaTermSamples, 1, quantile, probs = alpha / 2)
  gammaTermHi <- apply(gammaTermSamples, 1, quantile, probs = 1 - alpha / 2)

  dfList2[[k]] <- data.frame(
    j = terms[k],
    xj = x[, k],
    fjxj = gammaTermMean,
    lo = gammaTermLo,
    hi = gammaTermHi
  )
  
  print(k)
}

gammaDf2 <- dfList2 %>% plyr::rbind.fill()

## Plot it 
gammaDf2 %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_rug(aes(xj, fjxj), sides = "b", alpha = 0.1) +
  ## geom_line(aes(xj, lo), col = "firebrick3", lty = "dashed") +
  ## geom_line(aes(xj, hi), col = "firebrick3", lty = "dashed") +
  geom_ribbon(aes(xj, ymin = lo, ymax = hi), fill = "grey60", alpha = 0.5) +
  geom_line(aes(xj, fjxj), col = "firebrick3") +
  facet_wrap(~j, scales = "free_x") +
  labs(title = "Projected additive summary of GP fit",
       x = "",
       y = TeX("$f_j(x_j)$")) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

## Compare to previous fit
gamcompPlot <- rbind(gammaDf %>% mutate(mode = "without LAT/LON interaction"),
      gammaDf2 %>% mutate(mode = "with LAT/LON interaction")) %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_rug(aes(xj, fjxj), sides = "b", alpha = 0.1) +
  geom_line(aes(xj, lo, col = mode, group = mode), lty = "dashed", size = 1) +
  geom_line(aes(xj, hi, col = mode, group = mode), lty = "dashed", size = 1) +
  ## geom_ribbon(aes(xj, ymin = lo, ymax = hi), fill = "grey60", alpha = 0.5) +
  geom_line(aes(xj, fjxj, col = mode, group = mode), size = 1) +
  facet_wrap(~j, scales = "free_x") +
  scale_color_brewer("Summary", palette = "Dark2") + 
  labs(title = "Projected additive summaries of GP fit",
       subtitle = "Comparing summaries with and without LATITUDE/LONGITUDE interaction",
       x = "",
       y = TeX("$f_j(x_j)$")) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = "bottom") 

gamcompPlot

ggsave(here("figures/GP-gam-plot2.pdf"), gamcompPlot, width = 11.2, height = 8.8)

## Make the map

## If you want the whole rectangle...
## LATseq <- seq(min(califDf$LATITUDE), max(califDf$LATITUDE), length.out = 100)
## LONseq <- seq(min(califDf$LONGITUDE), min(califDf$LONGITUDE),  length.out = 100)

## mygrid <- expand.grid(LONseq, LATseq)

gamterms2 <- predict(gamFit2, newdata = califDf, type = "terms")

head(gamterms2)

## Plot

## Bounds of map
mapdata <- califDf %>%
  mutate(fLONLAT = gamterms2[, 4]) %>%
  rename(lon = LONGITUDE, lat = LATITUDE)

CA <- c(left = min(mapdata$lon),
        right = max(mapdata$lon),
        bottom = min(mapdata$lat),
        top = max(mapdata$lat))

library(ggmap)

map <- get_stamenmap(CA, zoom = 6, maptype = "toner-lite")

ggmap(map) +
  geom_point(data = mapdata, aes(lon, lat, col = fLONLAT), alpha = 0.3) +
  scale_color_viridis_c(TeX("$f(x_k, x_l)$"), option="A") + 
  labs(title = "Spatial trend",
       x = "", y = "")

ggsave(here("figures/GP-gam-plot2-spatial.pdf"), width = 7, height = 7)

###############################################################################
                                        #       GAM - yhat residuals (2)      #
###############################################################################

gamres2 <- as.numeric(yhat - gamma2)

resGAMtree2 <- rpart(gamres2 ~ ., data = califDf,
                    control = rpart.control(maxdepth=5, cp = 0.0001))

resGAMtree2 <- rpart(gamres2 ~ ., data = califDf)

rpart.plot(resGAMtree2,
           box.palette = "Greys", main = "Summarization residuals of additive summary fit to GP posterior mean\nafter adding LAT/LON interaction",
           cex = 1)

pdf(here("figures/GP-gam-res2.pdf"), width = 13, height = 5.7)
rpart.plot(resGAMtree2,
           box.palette = "Greys", main = "Residuals of GAM fit to GP posterior mean\nafter adding LAT/LON interaction",
           cex = 1)
dev.off()

###############################################################################
                                        #           LAT/LON, LAT/MHI          #
###############################################################################


gamFit311 <- gam(yhat ~
                 s(POPULATION) +
                 s(Median_rooms) +
                 s(LATITUDE, LONGITUDE) + 
                 s(LATITUDE, Median_household_income),
              data = califDf)

gamma311 <- predict(gamFit311, newdata = califDf, type = "response")

gamres311 <- as.numeric(y - gamma311)
resGAMtree311 <- rpart(gamres311 ~ ., data = califDf,
                    control = rpart.control(maxdepth=5, cp = 0.0001))

## rpart.plot(resGAMtree311,
##            box.palette = "Greys", main = "Residuals of GAM fit to GP posterior mean\nafter adding LAT/LON interaction",
##            cex = 1)


Xs311 <- model.matrix(gamFit311)

V311 <- vcov(gamFit311, dispersion = 1)
q311 <- crossprod(V311, crossprod(Xs311, yhat))
Q311 <- crossprod(V311, crossprod(Xs311, yhatSamples))

gammaSamples311 <- Xs311 %*% Q311

###############################################################################
                                        #           LAT/LON, LAT/MHI          #
###############################################################################


gamFit312 <- gam(yhat ~
                 s(POPULATION) +
                 s(Median_rooms) +
                 s(LATITUDE, LONGITUDE) + 
                 s(LONGITUDE, Median_household_income),
              data = califDf)

gamma312 <- predict(gamFit312, newdata = califDf, type = "response")

gamres312 <- as.numeric(y - gamma312)
resGAMtree312 <- rpart(gamres312 ~ ., data = califDf,
                    control = rpart.control(maxdepth=5, cp = 0.0001))


Xs312 <- model.matrix(gamFit312)

V312 <- vcov(gamFit312, dispersion = 1)
q312 <- crossprod(V312, crossprod(Xs312, yhat))
Q312 <- crossprod(V312, crossprod(Xs312, yhatSamples))

gammaSamples312 <- Xs312 %*% Q312

###############################################################################
                                        #           LAT/LON/MHI          #
###############################################################################


gamFit32 <- gam(yhat ~
                 s(POPULATION) +
                 s(Median_rooms) +
                 s(LONGITUDE, LATITUDE, Median_household_income),
              data = califDf)

gamma32 <- predict(gamFit32, newdata = califDf, type = "response")

gamres32 <- as.numeric(y - gamma32)
resGAMtree32 <- rpart(gamres32 ~ ., data = califDf,
                    control = rpart.control(maxdepth=5, cp = 0.0001))


Xs32 <- model.matrix(gamFit32)

V32 <- vcov(gamFit32, dispersion = 1)
q32 <- crossprod(V32, crossprod(Xs32, yhat))
Q32 <- crossprod(V32, crossprod(Xs32, yhatSamples))

gammaSamples32 <- Xs32 %*% Q32



###############################################################################
                                        #     Add interaction with MHI...     #
###############################################################################

## Three-way interaction...

gamFit31 <- gam(yhat ~
                 s(POPULATION) +
                 s(Median_rooms) +
                 s(LATITUDE, LONGITUDE, Median_household_income),
              data = califDf)

gamma31 <- predict(gamFit31, newdata = califDf, type = "response")

gamres31 <- as.numeric(y - gamma31)
resGAMtree31 <- rpart(gamres31 ~ ., data = califDf,
                    control = rpart.control(maxdepth=5, cp = 0.0001))


Xs31 <- model.matrix(gamFit31)

V31 <- vcov(gamFit31, dispersion = 1)
q31 <- crossprod(V31, crossprod(Xs31, yhat))
Q31 <- crossprod(V31, crossprod(Xs31, yhatSamples))

gammaSamples31 <- Xs31 %*% Q31


rpart.plot(resGAMtree31,
           box.palette = "Greys", main = "Residuals of GAM fit to GP posterior mean\nafter adding LAT/LON interaction",
           cex = 1)

pdf(here("GP/img2/GP-gam-res31.pdf"), width = 12, height = 12)
rpart.plot(resGAMtree31,
           box.palette = "Greys", main = "Residuals of GAM fit to GP posterior mean\nafter adding LAT/LON/MHI interaction",
           cex = 1)
dev.off()

###############################################################################
                                        #         Other interactions..        #
###############################################################################



mynames <- c("GAM with LAT,Income itx",
             "GAM with LAT,Rooms itx",
             "GAM with LON,Income itx",
             "GAM with LON,Rooms itx",
             "GAM with Income,Rooms itx")

## (3) LAT, MHI
gamFit3 <- gam(yhat ~
                 s(POPULATION) +
                 s(Median_rooms) +
                 s(LATITUDE, Median_household_income) +
                 s(LONGITUDE),
              data = califDf)

gamma3 <- predict(gamFit3, newdata = califDf, type = "response")

Xs3 <- model.matrix(gamFit3)

V3 <- vcov(gamFit3, dispersion = 1)
q3 <- crossprod(V3, crossprod(Xs3, yhat))
Q3 <- crossprod(V3, crossprod(Xs3, yhatSamples))

gammaSamples3 <- Xs3 %*% Q3

## (4) LAT, Med
gamFit4 <- gam(yhat ~
                s(POPULATION) +
                s(Median_household_income) +
                s(LONGITUDE) +
		s(LATITUDE, Median_rooms),
              data = califDf)

gamma4 <- predict(gamFit4, newdata = califDf, type = "response")

Xs4 <- model.matrix(gamFit4)

V4 <- vcov(gamFit4, dispersion = 1)
q4 <- crossprod(V4, crossprod(Xs4, yhat))
Q4 <- crossprod(V4, crossprod(Xs4, yhatSamples))

gammaSamples4 <- Xs4 %*% Q4


## (5) LON, MHI
gamFit5 <- gam(yhat ~
                s(POPULATION) +
                s(Median_rooms) +
                s(LONGITUDE, Median_household_income) +
		s(LATITUDE),
              data = califDf)

gamma5 <- predict(gamFit5, newdata = califDf, type = "response")

Xs5 <- model.matrix(gamFit5)

V5 <- vcov(gamFit5, dispersion = 1)
q5 <- crossprod(V5, crossprod(Xs5, yhat))
Q5 <- crossprod(V5, crossprod(Xs5, yhatSamples))

gammaSamples5 <- Xs5 %*% Q5


## (6) LON, Med
gamFit6 <- gam(yhat ~
                s(POPULATION) +
                s(Median_household_income) +
                s(LONGITUDE, Median_rooms) +
		s(LATITUDE),
              data = califDf)

gamma6 <- predict(gamFit6, newdata = califDf, type = "response")

Xs6 <- model.matrix(gamFit6)

V6 <- vcov(gamFit6, dispersion = 1)
q6 <- crossprod(V6, crossprod(Xs6, yhat))
Q6 <- crossprod(V6, crossprod(Xs6, yhatSamples))

gammaSamples6 <- Xs6 %*% Q6


## (7) MHI, Med
gamFit7 <- gam(yhat ~
                s(POPULATION) +
                s(Median_household_income, Median_rooms) +
		s(LONGITUDE) +  
		s(LATITUDE),
              data = califDf)

gamma7 <- predict(gamFit7, newdata = califDf, type = "response")

Xs7 <- model.matrix(gamFit7)

V7 <- vcov(gamFit7, dispersion = 1)
q7 <- crossprod(V7, crossprod(Xs7, yhat))
Q7 <- crossprod(V7, crossprod(Xs7, yhatSamples))

gammaSamples7 <- Xs7 %*% Q7


###############################################################################
                                        #         Compute diagnostics         #
###############################################################################


sigmaSamples <- sqrt(sigma2Samples)

rsqLin <- rsqGamma(gammaLin, yhatSamples)
rsqGAM2 <- rsqGamma(gammaSamples2, yhatSamples)
rsqGAM3 <- rsqGamma(gammaSamples3, yhatSamples)
rsqGAM4 <- rsqGamma(gammaSamples4, yhatSamples)
rsqGAM5 <- rsqGamma(gammaSamples5, yhatSamples)
rsqGAM6 <- rsqGamma(gammaSamples6, yhatSamples)
rsqGAM7 <- rsqGamma(gammaSamples7, yhatSamples)

rsqGAM311 <- rsqGamma(gammaSamples311, yhatSamples)
rsqGAM312 <- rsqGamma(gammaSamples312, yhatSamples)
rsqGAM32 <- rsqGamma(gammaSamples32, yhatSamples)


phiLin <- phiGamma(y, gammaLin, sqrt(sigma2Samples))
phiGAM3 <- phiGamma(y, gammaSamples3, sigmaSamples)
phiGAM4 <- phiGamma(y, gammaSamples4, sigmaSamples)
phiGAM5 <- phiGamma(y, gammaSamples5, sigmaSamples)
phiGAM6 <- phiGamma(y, gammaSamples6, sigmaSamples)
phiGAM7 <- phiGamma(y, gammaSamples7, sigmaSamples)

phiGAM311 <- phiGamma(y, gammaSamples311, sqrt(sigma2Samples))
phiGAM312 <- phiGamma(y, gammaSamples312, sqrt(sigma2Samples))
phiGAM32 <- phiGamma(y, gammaSamples32, sqrt(sigma2Samples))




###############################################################################
                                        #       All summary diagnostics       #
###############################################################################

library(RColorBrewer)

display.brewer.pal(3, "Greens")
display.brewer.pal(8, "Blues")

mygreens <- brewer.pal(3, "Greens")
myblues <- brewer.pal(6, "Blues")

## R-sq
mynames <- c("LAT/Income",
             "LAT/Rooms",
             "LON/Income",
             "LON/Rooms",
             "Income/Rooms")


rsqDfBig <- rbind(
  data.frame(rsq = rsqLin, model = "Linear", Name = "(NA)"),
  data.frame(rsq = rsqGAM, model = "Additive, with\nno itx", Name = "(NA)"),
  ## One 2-way itx
  data.frame(rsq = rsqGAM2, model = "Additive, with\none 2-way itx", Name = "LAT/LON"),
  data.frame(rsq = rsqGAM3, model = "Additive, with\none 2-way itx", Name = mynames[1]),
  data.frame(rsq = rsqGAM4, model = "Additive, with\none 2-way itx", Name = mynames[2]),
  data.frame(rsq = rsqGAM5, model = "Additive, with\none 2-way itx", Name = mynames[3]),
  data.frame(rsq = rsqGAM6, model = "Additive, with\none 2-way itx", Name = mynames[4]),
  data.frame(rsq = rsqGAM7, model = "Additive, with\none 2-way itx", Name = mynames[5]),
  ## Two 2-way itx
  data.frame(rsq = rsqGAM311,
             model = "Additive, with\ntwo 2-way itx",
             Name = "LAT/LON, LAT/Income"),
  data.frame(rsq = rsqGAM312,
             model = "Additive, with\ntwo 2-way itx",
             Name = "LAT/LON, LON/Income"),
  ## One 3-way interaction
  data.frame(rsq = rsqGAM32,
             model = "Additive, with\none 3-way itx",
             Name = "LAT/LON/Income")
)

theme_set(theme_bw(base_size = 18))

rsqPlot <- rsqDfBig %>%
  ggplot() +
  geom_violin(aes(model, rsq, col = Name), position = "identity", alpha = 0.1,
              draw_quantiles = 0.5, size = 0.75) +
  scale_color_manual("interactions\n(\"itx\")",
                     breaks = c("LAT/LON", mynames,
                                "LAT/LON, LAT/Income",
                                "LAT/LON, LON/Income",
                                "LAT/LON/Income"),
                     values = c("black",      #(NA)
                                "darkorange", #LAT/LON
                                myblues[2:6], #All other one 2-way
                                mygreens[2],  #GAM w/ LAT/LON,\nLAT/Income itx
                                mygreens[3],           #GAM w/ LAT/LON,\nLON/Income itx
                                "grey60")) +             #Three-way
  labs(title = "Variance of model explained by summarization",
       x = "Summary class", y = TeX("$R^2_\\gamma$")) +
  ## theme(legend.position = "bottom",
  ##       axis.title.y = element_text(vjust = 0.5, angle = 0, hjust = 0)) +
  theme(legend.position = "bottom"## ,
        ## axis.text.x = element_text(angle = 45, hjust = 1)
        )

rsqPlot

## phi 
phiDfBig <- rbind(
  data.frame(phi = phiLin, model = "Linear", Name = "(NA)"),
  data.frame(phi = phiGAM, model = "Additive, with\nno itx", Name = "(NA)"),
  ## One 2-way interaction
  data.frame(phi = phiGAM2, model = "Additive, with\none 2-way itx", Name = "LAT/LON"),
  data.frame(phi = phiGAM3, model = "Additive, with\none 2-way itx", Name = mynames[1]),
  data.frame(phi = phiGAM4, model = "Additive, with\none 2-way itx", Name = mynames[2]),
  data.frame(phi = phiGAM5, model = "Additive, with\none 2-way itx", Name = mynames[3]),
  data.frame(phi = phiGAM6, model = "Additive, with\none 2-way itx", Name = mynames[4]),
  data.frame(phi = phiGAM7, model = "Additive, with\none 2-way itx", Name = mynames[5]),
  ## Two 2-way interactions
  data.frame(phi = phiGAM311,
             model = "Additive, with\ntwo 2-way itx",
             Name = "LAT/LON, LAT/Income"),
  data.frame(phi = phiGAM312,
             model = "Additive, with\ntwo 2-way itx",
             Name = "LAT/LON, LON/Income"),
  ## One 3-way interaction
  data.frame(phi = phiGAM32,
             model = "Additive, with\none 3-way itx",
             Name = "LAT/LON/Income")
)

phiPlot <- phiDfBig %>%
  ggplot() +
  geom_violin(aes(model, phi, col = Name), position = "identity", alpha = 0.1,
              draw_quantiles = 0.5, size = 0.75) +
  scale_color_manual("interactions\n(\"itx\"",
                     breaks = c("LAT/LON", mynames,
                                "LAT/LON, LAT/Income",
                                "LAT/LON, LON/Income",
                                "LAT/LON/Income"),
                     values = c("black",      #(NA)
                                "darkorange", #LAT/LON
                                myblues[2:6], #All other one 2-way
                                mygreens[2],  #GAM w/ LAT/LON,\nLAT/Income itx
                                mygreens[3],           #GAM w/ LAT/LON,\nLON/Income itx
                                "grey60")) +             #Three-way
  labs(title = "Inflation of residual stdev after summarization",
       x = "Summary class", y = TeX("$\\phi_\\gamma$ (inverted scale)")) +
  ## theme(legend.position = "none", 
  ##       axis.title.y = element_text(vjust = 0.5, angle = 0)) +
  theme(legend.position = "none") + 
  scale_y_reverse()

phiPlot

                                        # Combine #############################
## Combine them
library(gridExtra)

mygrob <- arrangeGrob(rbind(ggplotGrob(rsqPlot), ggplotGrob(phiPlot), size="last"))

grid.arrange(mygrob)

ggsave("figures/GP-diagnostics-combined.pdf", plot = mygrob,
       width = 11.7, height = 12.8, units = "in")


