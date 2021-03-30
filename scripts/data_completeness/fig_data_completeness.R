## Generate Figure 1 of the paper: data completeness.
## By MW, GPLv3+, Oct. 2016
## The panels are the following:
## - i.   Map of the annual incidence in 2006
## - ii.  Map of the annual incidence in 2009
## - iii. Map of the annual incidence in 2012
## - iv.  Evolution of the number of HC and the incidence over time at the multi-district level
## - v.   Distribution of population mismatches
## - vi.  Distribution of cases mismatches
## - vii. Number of data points per HC
## - viii.Distribution of HC characteristic sizes over time

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

## ==== Reference stuffs ====
ref.pop <- read.csv("../data/brut/annuel-comparaison/Annuel_synthese.csv", sep=";")## Load the reference dataset to get the population references
ref.pop$district <- as.factor(tolower(as.character(ref.pop$district)))
ref.ld <- tolower(c("BOULSA", "DAFRA", "DANDE", "DÉDOUGOU", "DÔ", "GOURCY", "HOUNDÉ", "KARANGASSO VIGUÉ", "LÉNA", "OUAHIGOUYA", "ORODARA", "SÉGUÉNÉGA", "TITAO", "YAKO"))
ref.pop <- subset(ref.pop, district %in% ref.ld) # remove the unused data
ref.pop$district <- droplevels(ref.pop$district) # Remove the unused levels
## Rename the levels so that we match our BDD
levels(ref.pop$district)[levels(ref.pop$district)=="dédougou"] <- "dedougou"
levels(ref.pop$district)[levels(ref.pop$district)=="dô"] <- "do"
levels(ref.pop$district)[levels(ref.pop$district)=="houndé"] <- "hounde"
levels(ref.pop$district)[levels(ref.pop$district)=="karangasso vigué"] <- "karangasso_vigue"
levels(ref.pop$district)[levels(ref.pop$district)=="léna"] <- "lena"
levels(ref.pop$district)[levels(ref.pop$district)=="séguénéga"] <- "seguenega"

comp <- read.csv("../data/brut/annuel-comparaison/weekly-extraction.csv") # Reference cases
comp.w <- paste("W", str_pad(1:53, 2, pad="0"), sep="")
comp.name <- c("year", "district", comp.w)
names(comp) <- comp.name
comp.years <- 2004:2014
comp <- subset(comp, year %in% comp.years)
comp <- subset(comp, tolower(district) %in% c(ref.ld, "lena", "seguenega"))
comp.m <- melt(comp, id.vars=c("district", "year"))
names(comp.m) <- c("district", "year", "week", "cas")
comp.m$datew <- paste(comp.m$year, comp.m$week, "7", sep="-")
comp.m$date <- ISOweek2date(comp.m$datew)
comp.o <- melt(as.data.frame(menin.cas))
names(comp.o) <- c("nomfus", "date", "cas")
comp.o$district <- unlist(lapply(strsplit(comp.o$nomfus, "__"), function(i){i[[1]]}))
comp.m$district <- droplevels(factor(tolower(comp.m$district)))
levels(comp.m$district)[levels(comp.m$district)=="dédougou"] <- "dedougou"
levels(comp.m$district)[levels(comp.m$district)=="dô"] <- "do"
levels(comp.m$district)[levels(comp.m$district)=="houndé"] <- "hounde"
levels(comp.m$district)[levels(comp.m$district)=="karangasso vigué"] <- "karangasso_vigue"
levels(comp.m$district)[levels(comp.m$district)=="séguénéga"] <- "seguenega"

## ==== Functions
map.cols <- function(v, max.inc=NULL, colormap=viridis, reverse=FALSE, overflow="#ff0000") {
    if (is.null(max.inc)) {
        max.inc <- as.integer(max(v, na.rm=TRUE))
    }
    v.int <- as.integer(v)
    cb <- colormap(max.inc+1)
    if (reverse) {
        cb <- rev(cb)
    }
    cb <- c(cb, overflow)
    v.int[v.int>max.inc] <- max.inc+2
    cs <- cb[v.int+1]
    cs[is.na(cs)] <- "gray"
    return(list(col=cs,
                cb=cb,
                min=0,
                max=max.inc))
    }
}

# Function to plot color bar (from: http://www.colbyimaging.com/wiki/statistics/color-bars)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', new=TRUE)
{
    scale = (length(lut)-1)/(max-min)
    if (new) {
        dev.new(width=1.75, height=5)
    }
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

## 
## ==== Panel i-iii (new version) ====
##
dat <- c("2006-03-05", "2006-03-12", "2006-03-19", "2006-03-26", "2006-04-02")
inc2006 <- list()

for (d in dat) {
    inc2006[[d]] <- seasons.inc[["2006"]][,c("nomfus", d)]
    names(inc2006[[d]]) <- c("NAME", "inc")
    inc2006[[d]] <- merge(partition.merg[["2006"]], inc2006[[d]])
}


##par(mfrow=c(2,2))
for (d in dat[1:5]) {
    pdf(paste(fig.path, paste("weekly-incidence-", d, ".pdf", sep=""), sep="/"), width=12, height=10)
    colbar <- map.cols(inc2006[[d]]$inc*1e5, colormap=heat.colors, reverse=T, max.inc=200) # Get the colorbar object
    plot(bf.districts, main=d)
    plot(inc2006[[d]], col=colbar[["col"]], border="black", add=T) # Plot the data
    plot(bf.rivers, add=T, col=alpha("#0093AF", 0.3))
    plot(bf.roads, add=T, col=alpha("#FF9966", 0.6), lwd=2)
    plot(bf.country, lwd=2, add=T)
    plot(bf.cities, add=T, pch=20, cex=2, col="#ED2939")
    dev.off()
}

pdf(paste(fig.path, "weekly-incidence-colorbar.pdf", sep="/"), width=1.75, height=5)
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=FALSE)      # And the colorbar
dev.off()

##
## ==== Fig 3:
##
dat <- "2008-01-20"
inc2008 <- seasons.inc[["2008"]][,c("nomfus", dat)]
names(inc2008) <- c("NAME", "inc")
inc2008 <- merge(partition.merg[["2008"]], inc2008)


colbar <- map.cols(inc2008$inc*1e5, colormap=heat.colors, reverse=T, max.inc=NULL)

pdf(paste(fig.path, "clusturing-example.pdf", sep="/"), width=12, height=10)
plot(bf.districts, main=dat)
plot(inc2008, col=colbar[["col"]], border="black", add=T) # Plot the data
plot(bf.rivers, add=T, col=alpha("#0093AF", 0.3))
plot(bf.roads, add=T, col=alpha("#FF9966", 0.6), lwd=2)
plot(bf.country, lwd=2, add=T)
plot(bf.cities, add=T, pch=20, cex=2, col="#ED2939")
dev.off()

pdf(paste(fig.path, "clusturing-example-colorbar.pdf", sep="/"), width=1.75, height=5)
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=T)      # And the colorbar
dev.off()
## 
## ==== Panel i-iii (old version) ====
##

yearly.inc <- lapply(seasons.inc, function(l){
    l$inc <- rowSums(l[,1:(ncol(l)-2)], na.rm=TRUE);
    out <- l[,c("inc", "nomfus")]
    names(out) <- c("inc", "NAME") # To match the merge requirements
    return(out);
}) # Compute yearly incidences

yearly.shp <- lapply(names(yearly.inc), function(y) {
    return(merge(partition.merg[[y]],yearly.inc[[y]]))
}) # Merge the spatial and epidemiological data
names(yearly.shp) <- names(yearly.inc)

par(mfrow=c(2,2))
for (y in c(2006, 2009, 2012))
{
    pdf(paste(fig.path, paste("annual-incidence-", y, ".pdf", sep=""), sep="/"), width=12, height=10)
    y <- as.character(y)
    colbar <- map.cols(yearly.shp[[y]]$inc*1e5, colormap=viridis, reverse=F, max.inc=400) # Get the colorbar object
    
    plot(bf.districts, main=y)
    plot(yearly.shp[[y]], col=colbar[["col"]], border="black", add=T) # Plot the data
    plot(bf.rivers, add=T, col=alpha("#0093AF", 0.3))
    plot(bf.roads, add=T, col=alpha("#FF9966", 0.6), lwd=2)
    plot(bf.country, lwd=2, add=T)
    plot(bf.cities, add=T, pch=20, cex=2, col="#ED2939")
    dev.off()
}
pdf(paste(fig.path, "annual-incidence-colorbar.pdf", sep="/"), width=1.75, height=5)
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=FALSE)      # And the colorbar
dev.off()

## 
## ==== Panel iv ====
##

## Compute the data
agg.dat <- names(menin.cas[,-ncol(menin.cas)])
agg.cas <- colSums(menin.cas[,-ncol(menin.cas)], na.rm=T)
agg.pop <- colSums(menin.pop[,-ncol(menin.pop)], na.rm=T)
agg.tot <- data.frame(dat=as.Date(agg.dat), cas=agg.cas, pop=agg.pop, inc=agg.cas/agg.pop)

fst.dat <- as.Date(paste(names(seasons.inc)[1:(length(seasons.inc)-1)], "01", "01", sep="-"))
fst.count <- unlist(lapply(seasons.inc[1:(length(seasons.inc)-1)], function(df) {
    sum(as.numeric(unlist(apply(df[,1:(ncol(df)-2)],1, function(l) {!all(is.na(l))}))), na.rm=T)
}))
fst.tot <- data.frame(dat=fst.dat, fs=fst.count)

## Plot
pdf(paste(fig.path, "incidence-number-of-hc.pdf", sep="/"))
par(mar=c(5,4,4,5)+.1)
plot(fst.tot, ylim=c(0,max(fst.tot$fs)+30), pch=17, col="#FF9966", cex=2, ylab="Number of health centers", xlab="time")
lines(fst.tot, col="#FF9966", lwd=2)

par(new=TRUE)
plot(inc*1e5~dat, data=agg.tot,type="l",lwd=3, col="#1034A6",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext(expression(paste(plain("cases/10") ^5, plain(" inhabitants"), sep="")),side=4,line=3)

legend("topright", legend=c("Weekly incidence", "Number of health centers"),
       lwd=c(3,2),
       col=c("#1034A6", "#FF9966"),
       bty="n",
       )
dev.off()

##
## ==== Panel v ====
##

## ## Aggregate the data the way we want
## ref.agg <- lapply(seasons.pop, function(s) {
##     da <- s[,c(10, ncol(s)-1)]
##     names(da) <- c("pop", "nomfus")
##     da$district <- as.character(unlist(lapply(strsplit(da$nomfus, "__"), function(e){return(e[[1]])})))
##     return(aggregate(da$pop, by=list(da$district), FUN=function(ll){sum(ll, na.rm=T)}))
## })# Compute the population per district and season by summing the populations of the HC
## ref.agg <- lapply(names(ref.agg), function(y) {
##     names(ref.agg[[as.character(y)]]) <- c("district", "pop")
##     ref.agg[[y]]$season <- y
##     return(ref.agg[[y]])
## })
## names(ref.agg) <- names(seasons.pop)
## ref.bind <- do.call("rbind", ref.agg)
## ref.bind <- subset(ref.bind, pop!=0)
## ref.bind$district <- as.factor(ref.bind$district)
## names(ref.pop) <- c("district", "ref.pop", "season")

## ## all(as.character(levels(ref.pop$district))==as.character(levels(ref.bind$district))) # sanity chk
## ref.comp <- merge(ref.pop, ref.bind, by=c("district", "season"))
## ref.comp$err <- (ref.comp$pop-ref.comp$ref.pop)/ref.comp$ref.pop

## ## Plot
## ret <- barplot(sort(ref.comp$err), col='#545AA7')
## thr <- .15 
## lines(c(0, ret[length(ret)]), c(thr,thr), col="red", lty=2, lwd=3)
## lines(c(0, ret[length(ret)]), c(-thr,-thr), col="red", lty=2, lwd=3)

## ## Se what's wrong
## ref.comp[abs(ref.comp$err) > thr,]

## Other method
a <- read.csv("../data/brut/150317_Collecte_Districts/Annuel_synthese.csv", sep=';')
a$district <- as.factor(toupper(as.character(a$district)))
a$population <- as.numeric(as.character(a$population))
r <- read.csv("../data/brut/annuel-comparaison/Annuel_synthese.csv", sep=";")
r$district <- as.factor(toupper(as.character(r$district)))
ld <- c("BOULSA", "DAFRA", "DANDE", "DÉDOUGOU", "DÔ", "GOURCY", "HOUNDÉ", "KARANGASSO VIGUÉ", "LÉNA", "OUAHIGOUYA", "ORODARA", "SÉGUÉNÉGA", "TITAO", "YAKO")
rr <- subset(r, district %in% ld)
rr$district <- droplevels(rr$district)
ag <- aggregate(a$population, by=list(a$district, a$annee), FUN=function(p){sum(p, na.rm=T)})
names(ag) <- c("district", "annee", "population")
cmp <- merge(ag, r, by=c("district", "annee"), all.x=T)
dif <- cmp$population.x-cmp$population.y
er <- dif/cmp$population.x
res <- cbind(cmp,dif,er)
res <- res[!is.na(res$er),]

## Plot
pdf(paste(fig.path, "error-populations.pdf", sep="/"), width=10, height=4)
thr <- .025
ret <- barplot(sort(res$er), col='#545AA7',
               ylim=c(-.1, .1), xlab="Health-center-year", ylab="Relative error")
lines(c(0, ret[length(ret)]), c(thr,thr), col="red", lty=2, lwd=2)
lines(c(0, ret[length(ret)]), c(-thr,-thr), col="red", lty=2, lwd=2)
dev.off()

sum(abs(res$er) < thr)/length(res$er)
sum(abs(res$er) < .05)/length(res$er)
## Se what's wrong
res[abs(res$er) > thr,]

##
## ==== Panel vi ====
##
## tolower(unique(comp$district)) %in% as.character(comp.o$district)
comp.agg <- aggregate(comp.o$cas, by=list(comp.o$district, comp.o$date), FUN=function(i){sum(i, na.rm=TRUE)})
comp.1 <- aggregate(comp.o$cas, by=list(comp.o$district, comp.o$date), FUN=function(i){
    if (all(is.na(i))) {
        return(NA)
    } else {
        return(1)
    }})
comp.agg$x <- comp.agg$x * comp.1$x
names(comp.agg) <- c("district", "date", "cas.ref")
comp.agg$year <- substr(comp.agg$date, 1, 4)
comp.m$cas <- as.numeric(as.character(comp.m$cas))
comp.m.agg <- aggregate(comp.m$cas, by=list(comp.m$district, comp.m$year), FUN=function(i){sum(i, na.rm=TRUE)})
comp.m.1 <- aggregate(comp.m$cas, by=list(comp.m$district, comp.m$year), FUN=function(i){
    if (all(is.na(i))) {
        return(NA)
    } else {
        return(1)
    }})
comp.m.agg$x <- comp.m.agg$x * comp.m.1$x
names(comp.m.agg) <- c("district", "year", "cas")
comp.merg <- merge(comp.m, comp.agg, by=c("district", "date")) # agg is the reference
comp.merg$cas.ref <- as.numeric(as.character(comp.merg$cas.ref))
comp.merg$cas <- as.numeric(as.character(comp.merg$cas))
comp.merg$diff <- comp.merg$cas.ref-comp.merg$cas

pdf(paste(fig.path, "error-cases.pdf", sep="/"), width=10, height=4)
thr <- 2
ret <- barplot(sort(comp.merg$diff), col='#545AA7',
               ylim=c(-40, 40),
               xlab="Health-center-week",
               ylab="Absolute error (cases)",
               border=NA)
lines(c(0, ret[length(ret)]), c(thr,thr), col="red", lty=2, lwd=2)
lines(c(0, ret[length(ret)]), c(-thr,-thr), col="red", lty=2, lwd=2)
dev.off()
median(abs(comp.merg$diff), na.rm=T)
sum(abs(comp.merg$diff)<=5, na.rm=T)/sum(!is.na(comp.merg$diff))
##
## ==== Panel vii ====
##
dtp.present <- rowSums(!is.na(menin.inc[,1:(ncol(menin.inc)-1)]))
dtp.tab <- data.frame(nfs=dtp.present, NAME=menin.inc$nomfus, nyp=dtp.present/52.)
dtp.shp <- merge(partition.merg[["2012"]],dtp.tab)

## Plot
pdf(paste(fig.path, "number-data-points.pdf", sep="/"), width=12, height=10)
colbar <- map.cols(dtp.shp$nyp, colormap=viridis, reverse=F) # Get the colorbar object
    
plot(bf.districts)
plot(dtp.shp, col=colbar[["col"]], border="black", add=T) # Plot the data
plot(bf.rivers, add=T, col=alpha("#0093AF", 0.3))
plot(bf.roads, add=T, col=alpha("#FF9966", 0.6), lwd=2)
plot(bf.country, lwd=2, add=T)
plot(bf.cities, add=T, pch=20, cex=2, col="#ED2939")

dev.off()
pdf(paste(fig.path, "number-data-points-colorbar.pdf", sep="/"), width=1.75, height=5)

color.bar(colbar[["cb"]], colbar[["min"]], nticks=12, colbar[["max"]], new=F)      # And the colorbar

dev.off()

##
## ==== Panel viii ====
##
hc2004 <- partition.merg[["2004"]] ## Subset the data
hc2012 <- partition.merg[["2012"]]
proj4string(hc2004) <- CRS("+init=epsg:4326") ## Load the right projection
proj4string(hc2012) <- CRS("+init=epsg:4326")
hc2004 <- spTransform(hc2004, CRS("+init=epsg:32630")) ## Convert to metric units
hc2012 <- spTransform(hc2012, CRS("+init=epsg:32630"))

hc2004.area <- sapply(slot(hc2004, "polygons"), slot, "area")/1e6 ## get km^2
hc2012.area <- sapply(slot(hc2012, "polygons"), slot, "area")/1e6 ## get km^2

h2004 <- hist(hc2004.area**.5, breaks=seq(0,30, by=.5), plot=F)
h2012 <- hist(hc2012.area**.5, breaks=seq(0,30, by=.5), plot=F)

plot(h2004$mids, h2004$density, type='l')
lines(h2004$mids, h2012$density, col='red')

pdf(paste(fig.path, "hc-size.pdf", sep="/"), width=6, height=4)
bw=.75
plot(density(hc2004.area**.5, bw=bw),
     xlab="Characteristic radius (km)",
     lwd=4, col='gray', main='')
lines(density(hc2012.area**.5, bw=bw), lwd=4, lty=2)
legend("topright", c("2004", "2012"), col=c('gray', 'black'), lty=c(1,2), lwd=4, bty="n")
dev.off()
