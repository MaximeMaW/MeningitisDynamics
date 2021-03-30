## Generates a few statistics about the database
## By MW, GPLv3+, Oct. 2016
## The computed statistics are the following

## ==== Introduction: loading stuff
library(scales) # For transparent colors
library(rgdal)
library(rgeos)
library(maptools)
library(zoo)
library(fields)
library(ISOweek)
library(ncdf4)
library(reshape2)
library(stpp)
library(viridis) # Colormap
library(parallel) ## Code can e easily refactord to work without
library(stringr)
source("utilitaires.R")


fig.path = "../papers/woringer2/figures/data_completeness/"
amh <- read.csv("../data/epi/150413_bdd_FRFE.csv")
amh$cas[amh$district=="titao" & amh$date=="2011-03-27"] <- 0 ## Correction d'un bug pesant
menin <- subset(amh, maladie=="menin")

fus <- read.csv("../data/geolocalisation/150413_fusion_fs_2500m.csv") ## Carte des fusions
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg
load("../data/foyers/160517_foyers.RData") # foyers, comp, compp
basename <- "../data/carto"
bf.maps <- readOGR(dsn = paste(basename, "WEB", sep="/"), "BUF") # Main administrative info
bf.maps <- gBuffer(bf.maps, width=0, byid=TRUE) # The shapefile is somehow broken. Repair it.
bf.country <- unionSpatialPolygons(bf.maps, bf.maps$ADM0) # Subset the BF
bf.districts <- unionSpatialPolygons(bf.maps, bf.maps$ADM2)
bf.rivers <- readOGR(dsn = paste(basename, "BF", sep="/"), "RIVERS") # Main administrative info
bf.rivers <- gIntersection(bf.rivers, bf.country) # Make sure that we stay within the borders
bf.roads <- readOGR(dsn = paste(basename, "BF", sep="/"), "roads")
bf.roads <- subset(bf.roads, RDLNTYPE==1) # Keep only main roads
bf.roads <- gIntersection(bf.roads, bf.country) # Make sure that we stay within the borders
bf.cities <- readOGR(dsn = paste(basename, "BF", sep="/"), "airp") # Airports -> big cities

menin.cas <- fusionner.fs(menin, fus, "cas", verbose=FALSE)
menin.pop <- fusionner.fs(menin, fus, "population", verbose=FALSE)
menin.inc <- menin.cas
idx <- 1:(ncol(menin.inc)-1)
menin.inc[,idx] <- menin.cas[,idx]/menin.pop[,idx]
seasons.inc <- to.seasons(menin.inc) ## split into seasons
seasons.pop <- to.seasons(menin.pop) ## split into seasons

## ==== Computing various statistics on the data

## number of HC-weeks
sum(!is.na(as.data.frame(menin.cas[,1:(ncol(menin.cas)-1)])), na.rm=T) # 129342 (eq. 2487 hc-y)

## number of cases
sum(as.data.frame(menin.cas[,1:(ncol(menin.cas)-1)]), na.rm=T) # 15344

## Area covered
hc2012 <- partition.merg[["2012"]]
proj4string(hc2012) <- CRS("+init=epsg:4326")
hc2012 <- spTransform(hc2012, CRS("+init=epsg:32630"))
hc2012.area <- sapply(slot(hc2012, "polygons"), slot, "area")/1e6 ## get km^2
sum(hc2012.area)/274200 ## Area from wikipedia
## Answer: 0.2015859

## Median area
median(hc2012.area**.5)

## Max distance between two HC
max(dist(coordinates(hc2012)))**.5

## Population covered
pop.covered <- sum(menin.pop[,"2012-01-08"], na.rm=T)## 3671591
pop.frac <- pop.covered/16590813 # from http://data.worldbank.org/country/burkina-faso
## Answer: 0.2213027

##
## Early cases stats
##

## Number of early cases
mm <- menin
mm$month <- as.numeric(substr(mm$date, 6,7))
mm.s <- subset(mm, month >9 & as.Date(date)<"2014-08-01")
sum(mm.s$cas, na.rm=TRUE) # 471
sum(subset(menin, as.Date(date)<"2014-08-01")$cas, na.rm=TRUE) # number of cases, 15413

## Number of early cases before an epidemic event
cc <- subset(climep, seu==TRUE)
sum(cc[!duplicated(cc[,c("year", "nomfus")]),]$precoces, na.rm=TRUE) # 53

##
## Effect of the aggregation on the number of HC
##
length(menin.cas$nomfus) ## 362
dstr <- unlist(lapply(strsplit(menin.cas$nomfus, "--"), length))
sum(dstr) ## 456
