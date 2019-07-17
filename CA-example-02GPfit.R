
library(here)
library(ggplot2)
library(dplyr)
library(GPFDA)

source(here("CA-example-01dataprep.R"))

system.time(GPfit <- gpr(x, y))

## Covariance matrix and eigendecomposition
K <- cov.pow.ex(GPfit$hyper, GPfit$train.x)
Ke <- eigen(K)

save(list = ls(), file = "R-output/CA-GPfit.Rdata")



