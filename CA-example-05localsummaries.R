 
library(here)
library(ggplot2)
library(dplyr)
library(mgcv)
library(rpart)
library(rpart.plot)
## library(ggmap)
library(latex2exp)
library(splines)
## library(laGP)
library(GPFDA)
library(Matrix)
library(mvtnorm)
library(stringr)
library(RColorBrewer)
library(scales)
theme_set(theme_bw(base_size = 16))

library(tidycensus)
library(sf)
library(tigris)
census_api_key("3deb7c3e77d1747cf53071c077e276d05aa31407", install = TRUE, overwrite = TRUE)
library(rmapzen)
mz_set_tile_host_nextzen(key = ("hxNDKuWbRgetjkLAf_7MUQ"))

## Read in functions
source(here("R/rsq.R"))
source(here("R/phiGamma.R"))
source(here("R/rsqGamma.R"))
source(here("R/makeNewPoint.R"))
source(here("R/predictGPposterior.R"))
source(here("R/localLinearSummaryGP.R"))
source(here("R/localLinearSummaryGPone.R"))
source(here("R/posteriorCI.R"))

load(here("R-output/CA-GPfit.Rdata"))
load(here("R-output/GPgibbs-samples.Rdata"))




## Function for obtaining vector tiles
get_vector_tiles <- function(bbox){
  mz_box=mz_rect(bbox$xmin,bbox$ymin,bbox$xmax,bbox$ymax)
  mz_vector_tiles(mz_box)
}

###############################################################################
                                        #          Preliminary stuff          #
###############################################################################


## Shape file for California
cashapeall <- get_acs(state = "CA", geography = "state", geometry = TRUE,
                      variables = "B19013_001")$geometry

## Shape files for all counties in California
cashape <- get_acs(state = "CA", geography = "county", geometry = TRUE,
                   variables = "B19013_001")

cashapetract <- get_acs(state = "CA", geography = "tract",
                        geometry = TRUE, variables = "B19013_001")

## Land areas for selected counties
sfshape <- (cashape %>% filter(NAME == "San Francisco County, California"))$geometry
smshape <- (cashape %>% filter(NAME == "San Mateo County, California"))$geometry

sfshape2 <- st_union(sfshape, smshape) #Combine SF and SM counties

lashape <- (cashape %>% filter(NAME == "Los Angeles County, California"))$geometry
ocshape <- (cashape %>% filter(NAME == "Orange County, California"))$geometry

## Combine LA area counties
lashape2 <- st_union(lashape, ocshape)

frshape <- (cashape %>% filter(NAME == "Fresno County, California"))$geometry

## Plot to make sure it works
testPlot1 <- ggplot() +
  geom_sf(data = sfshape2) +
  geom_sf(data = smshape, aes(fill = "San Matteo")) +
  geom_sf(data = sfshape, aes(fill = "San Francisco")) +
  scale_fill_brewer(palette = "Blues")

print(testPlot1)

library(RColorBrewer)

testPlot2 <- ggplot() +
  geom_sf(data = lashape2) +
  geom_sf(data = lashape, aes(fill = "LA")) +
  geom_sf(data = ocshape, aes(fill = "Orange")) +
  scale_fill_brewer("County", palette = "Oranges") +
  labs(title = "Counties in LA MSA")

print(testPlot2)
  

## load("R-output/county-shapes.Rdata")

###############################################################################
                                        #       Make a map of california      #
###############################################################################

cabbox <- st_bbox(cashapeall)

testPlotCA <- cashapeall %>%
  ggplot() +
  geom_sf() +
  coord_sf(xlim = c(cabbox$xmin, cabbox$xmax),
           ylim = c(cabbox$ymin, cabbox$ymax))

print(testPlotCA)

ca_vector_tiles <- get_vector_tiles(cabbox)
names(ca_vector_tiles)

ca_water <- as_sf(ca_vector_tiles$water)
ca_roads <- as_sf(ca_vector_tiles$roads)

ca_roads$kind %>% unique()

mycols1 <- c("palegreen4", "darkorange", "dodgerblue3")

metroMap <- cashapeall %>%
  ggplot() +
  geom_sf(data = ca_water %>%
            mutate(Area = st_area(geometry) %>% as.numeric()) %>%
            filter(Area > 1), #filter out points
          fill = "lightgrey") +
  geom_sf(data = ca_roads %>% filter(kind == "highway"), col = "grey30") +
  geom_sf(alpha = 0.1) +
  geom_sf(data = frshape, aes(fill = "Fresno"), alpha = 0.5) +
  geom_sf(data = sfshape, aes(fill = "SF+SM"), alpha = 0.5) +
  geom_sf(data = smshape, aes(fill = "SF+SM"), alpha = 0.5) +
  geom_sf(data = lashape, aes(fill = "LA+OC"), alpha = 0.5) +
  geom_sf(data = ocshape, aes(fill = "LA+OC"), alpha = 0.5) +
  coord_sf(xlim = c(cabbox$xmin, cabbox$xmax),
           ylim = c(cabbox$ymin, cabbox$ymax)) +
  theme_light(base_size = 16) +
  scale_fill_manual("Metro area", values = mycols1) + 
  labs(title = "Selected metropolitan areas in California") +
  theme(legend.position = "bottom")

## print(metroMap)

ggsave("figures/GP-CA-map.pdf", metroMap,
       width = 9, height = 8.5, units = "in")

###############################################################################
                                        #         Indices of counties         #
###############################################################################

smInd <- which(calif$County %in% c("San Mateo County"))
sfInd <- which(calif$County %in% c("San Francisco County"))

sfInd2 <- which(calif$County %in% c("San Francisco County", "San Mateo County"))

laInd <- which(calif$County == "Los Angeles County")
ocInd <- which(calif$County == "Orange County")

laInd2 <- which(calif$County %in% c("Los Angeles County", "Orange County"))

frInd <- which(calif$County == "Fresno County")

###############################################################################
                                        #          Dataframe subsets          #
###############################################################################

smdf <- calif %>% filter(County %in% c("San Mateo County"))
sfdf <- calif %>% filter(County %in% c("San Francisco County"))

sfdf2 <- calif %>% filter(County %in% c("San Francisco County", "San Mateo County"))

ladf <- calif %>% filter(County == "Los Angeles County")
ocdf <- calif %>% filter(County == "Orange County")

ladf2 <- calif %>% filter(County %in% c("Los Angeles County", "Orange County"))

frdf <- calif %>% filter(County == "Fresno County")

###############################################################################
                                        #              New points             #
###############################################################################

Nnew <- 1000

termsSub <- terms[1:3]

smSummary <- localLinearSummaryGP(smdf, smInd, smshape, Nnew, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "San Mateo")

sfSummary <- localLinearSummaryGP(sfdf, sfInd, sfshape, Nnew, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "San Francisco")

sf2Summary <- localLinearSummaryGP(sfdf2, sfInd2, sfshape2, Nnew, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "SF+SM")

laSummary <- localLinearSummaryGP(ladf, laInd, lashape, Nnew * 4, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "Los Angeles")

ocSummary <- localLinearSummaryGP(ocdf, ocInd, ocshape, Nnew, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "Orange County")

la2Summary <- localLinearSummaryGP(ladf2, laInd2, lashape2, Nnew * 4, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "LA+OC")

frSummary <- localLinearSummaryGP(frdf, frInd, frshape, Nnew, x, y,
                                  sigma2Samples, termsSub,
                                  Name = "Fresno County")

summaryDf <- rbind(
  smSummary$rsqGamma,
  sfSummary$rsqGamma,
  sf2Summary$rsqGamma,
  laSummary$rsqGamma,
  ocSummary$rsqGamma,
  la2Summary$rsqGamma,
  frSummary$rsqGamma
) 

summaryDf <- rbind(
  smSummary$rsqGamma %>% mutate(Level = "County"),
  sfSummary$rsqGamma %>% mutate(Level = "County"),
  sf2Summary$rsqGamma %>% mutate(Level = "MSA"),
  laSummary$rsqGamma %>% mutate(Level = "County"),
  ocSummary$rsqGamma %>% mutate(Level = "County"),
  la2Summary$rsqGamma %>% mutate(Level = "MSA"),
  frSummary$rsqGamma %>% mutate(Level = "County")
) 

## summaryDf %>%
##   ggplot() +
##   geom_density(aes(rsq, group = Name, col = Name), size = 1) +
##   scale_color_brewer(palette = "Dark2")

## summaryDf %>%
##   ggplot() +
##   geom_density(aes(rsq, group = Name, col = Name), size = 1) +
##   facet_wrap(~Level, ncol = 1) + 
##   scale_color_brewer(palette = "Dark2")

###############################################################################
                                        #  Try with just a few census tracts  #
###############################################################################

## Data for tracts for all of CA, then narrow to SF
cashapetract <- get_acs(state = "CA", geography = "tract",
                        geometry = TRUE, variables = "B19013_001", refresh=TRUE)

## cashapetract <- get_acs(state = "TX", geography = "block",
##                         geometry = TRUE, variables = "B19013_001")

## foobar <- read_sf("tl_2018_us_county.shp",)

## cashapetract <- get_decennial(state = "CA", geography = "tract",
##                         geometry = TRUE, variables = "P005003")

sfshapetract <- cashapetract %>%
  filter(str_detect(NAME, "San Francisco County"))

## Geographic features
bbox <- st_bbox(sfshapetract)
vector_tiles <- get_vector_tiles(bbox)
names(vector_tiles)

water <- as_sf(vector_tiles$water)
roads <- as_sf(vector_tiles$roads)

## Plot all tracts
SFtractPlot <- ggplot() +
  geom_sf(data = water, fill = "lightgrey") +
  geom_sf(data = roads, col = "grey30") +
  geom_sf(data=sfshapetract, aes(fill = estimate)) +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax),
           ylim = c(bbox$ymin, bbox$ymax)) +
  scale_fill_viridis_c()

SFtractPlot

## ## Not all tracts available
## sfsub <- sfshapetract %>%
##   filter(as.numeric(GEOID) %in% calif$GEO.id2)

## ggplot() +
##   geom_sf(data = water, fill = "lightgrey") +
##   geom_sf(data = roads %>% filter(kind %in% c("highway", "major_road")), col = "grey30") + 
##   geom_sf(data = sfsub, aes(fill = estimate)) + 
##   coord_sf(xlim = c(-122.55, bbox$xmax),
##            ylim = c(bbox$ymin, bbox$ymax)) +
##   scale_fill_viridis_c()

## mapview::mapview(sfsub)

###############################################################################
                                        #           Selected tracts           #
###############################################################################

## ## IDs for selected tracts
mygeoid1 <- c("06075017801", "06075017601", "06075017802",
              "06075018000", "06075017700", "06075022801",
              "06075022803", "06075020900", "06075022803",
              "06075022901", "06075022902", "06075022903")

mygeoid2 <- c("06075035201", "06075035202", "06075035100",
              "06075032700", "06075032602", "06075032601",
              "06075035400", "06075035300", "06075032902",
              "06075032802")

mygeoid3 <- c("06075047902", "06075047802", "06075047901",
              "06075047801", "06075042700", "06075047701",
              "06075047702", "06075045200")

as.numeric(mygeoid1) %in% calif$GEO.id2
as.numeric(mygeoid2) %in% calif$GEO.id2
as.numeric(mygeoid3) %in% calif$GEO.id2

## What the neighborhoods look like
sfsub1 <- sfshapetract %>%
  filter(GEOID %in% mygeoid1)

sfsub2 <- sfshapetract %>%
  filter(GEOID %in% mygeoid2)

sfsub3 <- sfshapetract %>%
  filter(GEOID %in% mygeoid3)

## Plot these locations
## ggplot() +
##   geom_sf(data = water, fill = "lightgrey") +
##   geom_sf(data = roads, col = "grey30") + 
##   geom_sf(data = sfsub1, aes(fill = "Area 1")) +
##   geom_sf(data = sfsub2, aes(fill = "Area 2")) +
##   geom_sf(data = sfsub3, aes(fill = "Area 3")) + 
##   coord_sf(xlim = c(-122.55, bbox$xmax),
##            ylim = c(bbox$ymin, bbox$ymax)) +
##   scale_fill_brewer("", palette = "Blues")

## Do indices, dataframe, summary...
sfsub1Ind <- which(calif$GEO.id2 %in% as.numeric(mygeoid1))
sfsub1Df <- calif %>% filter(GEO.id2 %in% as.numeric(mygeoid1))
sfsub1shape <- st_union(sfshapetract %>% filter(GEOID %in% mygeoid1))

sfsub2Ind <- which(calif$GEO.id2 %in% as.numeric(mygeoid2))
sfsub2Df <- calif %>% filter(GEO.id2 %in% as.numeric(mygeoid2))
sfsub2shape <- st_union(sfshapetract %>% filter(GEOID %in% mygeoid2))

sfsub3Ind <- which(calif$GEO.id2 %in% as.numeric(mygeoid3))
sfsub3Df <- calif %>% filter(GEO.id2 %in% as.numeric(mygeoid3))
sfsub3shape <- st_union(sfshapetract %>% filter(GEOID %in% mygeoid3))

sfsub1Summary <- localLinearSummaryGP(sfsub1Df, sfsub1Ind,
                                      sfsub1shape, Nnew, x, y,
                                      sigma2Samples, termsSub,
                                      Name = "SF Area 1")

sfsub2Summary <- localLinearSummaryGP(sfsub2Df, sfsub2Ind,
                                      sfsub2shape, Nnew, x, y,
                                      sigma2Samples, termsSub,
                                      Name = "SF Area 2")

sfsub3Summary <- localLinearSummaryGP(sfsub3Df, sfsub3Ind,
                                      sfsub3shape, Nnew, x, y,
                                      sigma2Samples, termsSub,
                                      Name = "SF Area 3")

sfsub1Summary$gamma %>% dim()

summaryDf <- rbind(
  smSummary$rsqGamma %>% mutate(Level = "County"),
  sfSummary$rsqGamma %>% mutate(Level = "County"),
  sf2Summary$rsqGamma %>% mutate(Level = "Metropolitan Area"),
  laSummary$rsqGamma %>% mutate(Level = "County"),
  ocSummary$rsqGamma %>% mutate(Level = "County"),
  la2Summary$rsqGamma %>% mutate(Level = "Metropolitan Area"),
  frSummary$rsqGamma %>% mutate(Level = "Metropolitan Area"),
  sfsub1Summary$rsqGamma %>% mutate(Level = "Neighborhood"),
  sfsub2Summary$rsqGamma %>% mutate(Level = "Neighborhood"),
  sfsub3Summary$rsqGamma %>% mutate(Level = "Neighborhood")
) %>%
  mutate(Level = factor(Level, levels = c("MSA", "County", "Neighborhood")))

## summaryDf %>%
##   ggplot() +
##   geom_density(aes(rsq, group = Name, col = Name), size = 1) +
##   facet_wrap(~Level, ncol = 1) +
##   labs(x = TeX("$\\R^{2}_{\\gamma}$")) +
##   geom_vline(xintercept = foobar)

###############################################################################
                                        #             Coefficients

                                        #
###############################################################################

## SOMA?
## https://data.sfgov.org/Geographic-Locations-and-Boundaries/Census-2010-Tracts-for-San-Francisco/rarb-5ahf/data

###############################################################################
                                        #          Just for one tract         #
###############################################################################

mytract <- mygeoid3[1]
mytractDf <- calif %>% filter(GEO.id2 == as.numeric(mytract))
mytractshape <- sfshapetract %>% filter(GEOID == mytract)

mytractsummary <- localLinearSummaryGPone(mytractDf, Nnew, x, y,
                                          sigma2Samples, termsSub, Name = "SF Tract")

summaryDf2 <- rbind(
  smSummary$rsqGamma %>% mutate(Level = "County"),
  sfSummary$rsqGamma %>% mutate(Level = "County"),
  sf2Summary$rsqGamma %>% mutate(Level = "Metropolitan Area"),
  laSummary$rsqGamma %>% mutate(Level = "County"),
  ocSummary$rsqGamma %>% mutate(Level = "County"),
  la2Summary$rsqGamma %>% mutate(Level = "Metropolitan Area"),
  frSummary$rsqGamma %>% mutate(Level = "Metropolitan Area"),
  sfsub1Summary$rsqGamma %>% mutate(Level = "Neighborhood"),
  sfsub2Summary$rsqGamma %>% mutate(Level = "Neighborhood"),
  sfsub3Summary$rsqGamma %>% mutate(Level = "Neighborhood"),
  mytractsummary$rsqGamma %>% mutate(Level = "Tract")
) %>%
  mutate(Level = factor(Level, levels = c("Metropolitan Area", "County", "Neighborhood", "Tract"))) %>%
  mutate(Name = factor(Name,
                       levels = c("LA+OC", "SF+SM",  "Fresno County",
                                  "Los Angeles", "Orange County" ,
                                  "San Mateo", "San Francisco",
                                  "SF Area 1", "SF Area 2",
                                  "SF Area 3", "SF Tract")))

blues <- brewer.pal(7, "Blues") %>% rev()
oranges <- brewer.pal(4, "Oranges") %>% rev()

mycols <- c(oranges[1], "navy", "palegreen4", oranges[2:3], blues[1:6])

diagPlot <- summaryDf2 %>%
  ggplot() +
  geom_density(aes(rsq, group = Name, col = Name), size = 1) +
  facet_wrap(~Level, ncol = 1) +
  labs(x = TeX("$\\R^{2}_{\\gamma}$")) +
  scale_color_manual(values = mycols)

diagPlot2x2 <- summaryDf2 %>%
  ggplot() +
  geom_density(aes(rsq, group = Name, col = Name), size = 1) +
  facet_wrap(~Level) +
  labs(x = TeX("$\\R^{2}_{\\gamma}$")) +
  scale_color_manual(values = mycols)

print(diagPlot)

print(diagPlot2x2)

ggsave("figures/GP-local-diagnostics-level.pdf", diagPlot,
       width = 7, height = 11, units = "in")

ggsave("figures/GP-local-diagnostics-level2x2.pdf", diagPlot2x2)

sfAreasMap <- ggplot() +
  geom_sf(data = water, fill = "lightgrey") +
  geom_sf(data = roads %>% filter(kind %in% c("highway", "major_road")), col = "grey30") + 
  geom_sf(data = sfsub1, aes(fill = "Area 1")) +
  geom_sf(data = sfsub2, aes(fill = "Area 2")) +
  geom_sf(data = sfsub3, aes(fill = "Area 3")) +
  geom_sf(data = mytractshape, aes(col = "Selected tract"), size = 1.5, alpha=0.1) +
  coord_sf(xlim = c(-122.55, bbox$xmax),
           ylim = c(bbox$ymin, bbox$ymax)) +
  ## scale_fill_brewer("", palette = "Set3") +
  scale_fill_brewer("", palette = "Blues") + 
  scale_color_manual("", values = "deeppink")

## print(sfAreasMap)

ggsave("figures/GP-SF-areas.png", sfAreasMap,
       width = 8.5, height = 7, units = "in")



###############################################################################
                                        #          Plot coefficients          #
###############################################################################

## For MSAs...
ciDf <- rbind(
  frSummary$ci %>% mutate(Level = "County", Name = "Fresno"),
  la2Summary$ci %>% mutate(Level = "MSA"),
  sf2Summary$ci %>% mutate(Level = "MSA")
  )

cityCoefPlot <- ciDf %>%
  ggplot() +
  geom_errorbar(aes(Name, ymin = lo, ymax = hi, col = Name), width = 0.5) +
  geom_point(aes(Name, est), col = "firebrick3") + 
  facet_wrap(~terms, scales = "free_x") +
  coord_flip() +
  scale_color_manual("",values=mycols1) +
  labs(title = "Local linear summaries of GP fit at metro area level",
       subtitle = "Between-city heterogeneity",
       x = "", y = "coefficient") +
  theme(legend.position = "none")

print(cityCoefPlot)

ggsave("figures/GP-city-coefs.pdf", cityCoefPlot,
       width = 9, height = 6.5, units = "in")

  ## Within SF...
ciDf2 <- rbind(
  smSummary$ci %>% mutate(Level = "County"),
  sfSummary$ci %>% mutate(Level = "County"),
  sf2Summary$ci %>% mutate(Level = "MSA"),
  sfsub1Summary$ci %>% mutate(Level = "Neighborhood"),
  sfsub2Summary$ci %>% mutate(Level = "Neighborhood"),
  sfsub3Summary$ci %>% mutate(Level = "Neighborhood"),
  mytractsummary$ci %>% mutate(Level = "Tract")
) %>%
  mutate(Level = factor(Level, levels = c("MSA", "County", "Neighborhood", "Tract"))) %>%
  mutate(Name = factor(Name,
                       levels = c("LA+OC", "SF+SM", "Los Angeles", "Orange County",
                                  "Fresno County", "San Mateo", "San Francisco",
                                  "SF Area 1", "SF Area 2", "SF Area 3", "SF Tract")))

SFcoefPlot <- ciDf2 %>%
  ggplot() +
  geom_errorbar(aes(Name, ymin = lo, ymax = hi, col = Level), width = 0.5) +
  geom_point(aes(Name, est), col = "firebrick3") + 
  facet_wrap(~terms, scales = "free_x") +
  coord_flip() +
  labs(title = "Local linear summaries of GP fit",
       subtitle = "Different levels of aggregation in SF",
       x = "", y = "coefficient") +
  scale_color_brewer(palette = "Dark2")

print(SFcoefPlot)

ggsave("figures/GP-SF-coefs.pdf", SFcoefPlot,
       width = 10.5, height = 8, units = "in")


