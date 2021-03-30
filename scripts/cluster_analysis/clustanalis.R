##
##      ==== Clustanalis ====
##
## A library to plot epidemic clusters in a "fancy" way
## By MW, Oct. 2016
## GPLv3+


## ==== Imports
library(rgdal)
library(scales) # For transparent colors
source("utilitaires.R")

library(maptools) # Only required for tests
library(rgeos) # idem

## ==== Functions
getSeason <- function(da){
    ## Return the season of a date
    y <- as.numeric(substr(da, 1, 4))
    if (da < as.Date(paste(y, 07, 01, sep="-"))) {
        return(y)
    } else {
        return(y+1)
    }
}

random.col <- function(nothing=NULL) {
    rgb(runif(1, .5), runif(1, .5), runif(1, .5), .3)
}

plot.cluster <- function(cls, cl.id, col="black",
                         dot.scale=1, plot.lines=TRUE, to.plot=c("time", "space"),
                         add=FALSE, fixaxis.y=NULL, fixaxis.t=NULL, thr.km=NULL,
                         sqrt.cas=TRUE) {
    ## Plots a simple cluster
    cl <- subset(cls, component==cl.id) ## Extract the cluster
    if (to.plot=="time"){
        cl$xx <- cl$x        
        cl$x <- cl$date
    }
        
    if (is.null(fixaxis.y)) {
        ylim=c(min(cl$y)-25e3, max(cl$y)+25e3)
    } else {
        ylim=c(min(cl$y)-.1*fixaxis.y, min(cl$y)+.9*fixaxis.y)
    }
    if (is.null(fixaxis.t)) {
        xlim=c(min(cl$date)-7, max(cl$date)+7)
    } else {
        xlim=c(min(cl$date)-.1*fixaxis.t, min(cl$date)+.9*fixaxis.t)
    }
        
    
    if (!add) {
        plot(y~x, data=cl, pch="", main=cl$season[1],
             ylim=ylim, xlim=xlim, xlab=to.plot)
    }
    if (plot.lines) {
        for (i in 1:nrow(cl)) {
            for (j in 1:nrow(cl)) {
                if (abs(cl[i,"date"]-cl[j,"date"])<=7) {
                    di <- ((cl[i,"xx"]-cl[j,"xx"])**2 +
                           (cl[i,"y"]-cl[j,"y"])**2)**0.5
                    if (is.null(thr.km) || di < thr.km*1000) {
                        lines(c(cl[i,"date"], cl[j,"date"]),
                              c(cl[i,"y"], cl[j,"y"]), col="lightgray")
                    }
                }
            }
        }
    }
    if (sqrt.cas) {
        cas.cex <- sqrt(cl$cas)*dot.scale
    } else {
        cas.cex <- cl$cas*dot.scale
    }
    points(y~x, data=cl,
         pch=16, cex=cas.cex, col=col)
}
stack.hist <- function(dat) {
    ## Generate the required format to plot a stack histogram
    ## Concretely, this guy makes a cumulative sum in the reverse
    ##+order so that (opaque) bars can overlap and that everything looks nice
    ## /!\ Plot with OPAQUE BARS

    tmp <- dat[,ncol(dat):1] ## Reverse input
    mask <- tmp/tmp
    
    res <- dat
    res[,] <- NA
    res[,1] <- tmp[,1]
    for (i in 2:ncol(dat)) {
        res[,i] <- rowSums(tmp[,1:i], na.rm=TRUE)*mask[,i]
    }
    return(res[,ncol(dat):1]) ## Reverse again
}

plot.timeline <- function(cas.ori, conn.all, conn.idx,
                          cols=rainbow(length(conn.idx)),
                          t.min=NULL, t.max=NULL) {
    ## Plots a simple timeline to locate the epidemic events wrt. the
    ##+general pattern of event
    
    ## ==== Parse arguments
    if (is.null(t.min)) {
        t.min <- min(as.Date(cas.all$date))
    }
    if (is.null(t.max)) {
        t.max <- max(as.Date(cas.all$date))
    }
    cas.all <- subset(cas.ori, date<t.max & date>=t.min)
    solid.cols <- substr(cols, 1, 7) ## /!\ careful if input is not in HEX format
    ## ==== Generate foreground params
    cas.new <- cas.all
    res.tmp <- cas.all
    res.tmp$cas <- NULL
    
    for (cl.id in conn.idx) {
        cl <- subset(conn.all, component==cl.id) ## Extract the cluster
        tl <- aggregate(cl$cas, by=list(cl$date), sum)
        names(tl) <- c("date", "cas")
        tl <- tl[order(tl$date),]
        tl$date <- as.Date(tl$date)
        cas.new[,paste("cluster", cl.id, sep="")] <- merge(res.tmp, tl, by="date", all.x=TRUE)$cas
    }
    ord <- order(colSums(cas.new[,3:ncol(cas.new)],na.rm=TRUE), decreasing=TRUE)
    cols <- cols[ord]
    solid.cols <- solid.cols[ord]
    sta <- stack.hist(cas.new[,2+ord])

    ## ==== Plot background
    with(cas.all,
         barplot(cas, names.arg=date, col="#DADADA", border=NA)
         )

    for (i in 1:(ncol(sta))) {
        barplot(sta[,i],
                col=solid.cols[i], border=NA, add=TRUE, ylim=c(0,10))
    }
    
}

## ==== Load data
load(file="../data/foyers/160526_foyers.RData") ## comp, compp, conn
fus <- read.csv("../data/geolocalisation/150413_fusion_fs_2500m.csv") ## Carte des fusions
amh <- read.csv("../data/epi/150413_bdd_FRFE.csv")
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg

basename <- "../data/carto"
bf.maps <- readOGR(dsn = paste(basename, "WEB", sep="/"), "BUF") # Main administrative info
bf.maps <- gBuffer(bf.maps, width=0, byid=TRUE) # The shapefile is somehow broken. Repair it.
bf.country <- unionSpatialPolygons(bf.maps, bf.maps$ADM0) # Subset the BF
bf.districts <- unionSpatialPolygons(bf.maps, bf.maps$ADM3)
bf.rivers <- readOGR(dsn = paste(basename, "BF", sep="/"), "RIVERS") # Main administrative info
bf.rivers <- gIntersection(bf.rivers, bf.country) # Make sure that we stay within the borders

proj4string(bf.country) <- CRS("+init=epsg:4326")
bf.country <- spTransform(bf.country, CRS("+init=epsg:32630")) # metric projection
proj4string(bf.districts) <- CRS("+init=epsg:4326")
bf.districts <- spTransform(bf.districts, CRS("+init=epsg:32630")) # metric projection
bf.rivers <- spTransform(bf.rivers, CRS("+init=epsg:32630")) # metric projection

## ==== Generate usual objects
amh$cas[amh$district=="titao" & amh$date=="2011-03-27"] <- 0 ## Correction d'un bug pesant
menin <- subset(amh, maladie=="menin")
menin.cas <- fusionner.fs(menin, fus, "cas", verbose=FALSE)
menin.pop <- fusionner.fs(menin, fus, "population", verbose=FALSE)
menin.inc <- menin.cas
idx <- 1:(ncol(menin.inc)-1)
menin.inc[,idx] <- menin.cas[,idx]/menin.pop[,idx]

##
## ==== Test code -- generate variables
##
cas.all <- data.frame(date=as.Date(names(menin.cas[1:(ncol(menin.cas)-1)])))
cas.all$cas <- colSums(menin.cas[,1:(ncol(menin.cas)-1)], na.rm=TRUE)

comp.list <- c(11, 13, 20, 31, 36, 62, 64, 81, 77, 82, 85, 96)
col.default <- rgb(0, 150, 175, 77, maxColorValue=255)
cols <- apply(data.frame(1:length(comp.list)), 1, FUN=random.col)
cols <- rainbow(length(comp.list), alpha=.3)
fixaxis.y <- 1.5e5
fixaxis.t <- 7*12

##
## ==== Test code -- timeline
##
pdf("../papers/woringer2/figures/cluster_analysis/timeline.pdf", width=14, height=4)
plot.timeline(cas.all, conn, comp.list,
              t.min="2006-01-01", t.max="2010-07-01")
dev.off()

##
## ==== Test code -- Plot temporal representation
##
pdf("../papers/woringer2/figures/cluster_analysis/clusters.pdf", width=14, height=8)

#plot.cluster(conn, comp.list[2], col=col.default, dot.scale=1, to.plot="time",
#             plot.lines=TRUE, add=FALSE, thr.km=10)

for (j in 1:length(comp.list)) {
    i <- comp.list[j]
    plot.cluster(conn, i, col=cols[j], dot.scale=1, to.plot="time",
                 add=FALSE, fixaxis.y=fixaxis.y, fixaxis.t=fixaxis.t, thr.km=14)
}

dev.off()

##
## ==== Test code -- Plot spatial representation
##

## Start with an ugly basemap
pdf("../papers/woringer2/figures/cluster_analysis/map.pdf", width=10, height=8)
plot(bf.districts, col="#FAFAFA")
plot(hc2012, add=T, col="white")

##plot.cluster(conn, comp.list[1], col=col.default, dot.scale=.25, to.plot="space",
##             plot.lines=FALSE, add=TRUE)

for (j in 1:length(comp.list)) {
    i <- comp.list[j]
    plot.cluster(conn, i, col=cols[j], dot.scale=.5, to.plot="space", 
                 add=TRUE, plot.lines=FALSE,
                 fixaxis.y=fixaxis.y, fixaxis.t=fixaxis.t)
}
dev.off()

## ==== Compute the areas of the shapefiles
partition.areas <- lapply(partition.merg, function(l) {
    proj4string(l) <- CRS("+init=epsg:4326")
    l.np <- spTransform(l, CRS("+init=epsg:32630")) # metric projection
    l.res <- data.frame(nomfus=as.character(l$NAME))
    l.res$area <- sapply(slot(l.np, "polygons"), slot, "area")/1e6 ## get km^2
    return(l.res)
})
partition.areas <- lapply(names(partition.areas), function(na){
    res <- partition.areas[[na]]
    res$season <- na
    return(res)
})
names(partition.areas) <- names(partition.merg)
partition.areas.df <- do.call("rbind", partition.areas)
partition.areas.df$nomfus <- as.character(partition.areas.df$nomfus)

## ==== Compute the centroids of the shapefiles
partition.centroids <- lapply(partition.merg, function(l) {
    proj4string(l) <- CRS("+init=epsg:4326")
    l.np <- spTransform(l, CRS("+init=epsg:32630")) # metric projection
    l.res <- data.frame(nomfus=as.character(l$NAME))
    co <- sapply(slot(l.np, "polygons"), function(ll){
        colMeans(slot(slot(ll, "Polygons")[[1]], "coords"))
    })
    l.res$x <- co[1,]
    l.res$y <- co[2,]
    return(l.res)
})

partition.centroids <- lapply(names(partition.centroids), function(na){
    res <- partition.centroids[[na]]
    res$season <- na
    return(res)
})
names(partition.centroids) <- names(partition.merg)
partition.centroids.df <- do.call("rbind", partition.centroids)
partition.centroids.df$nomfus <- as.character(partition.centroids.df$nomfus)


## ==== Generate the conn object (to be saved with comp, compp
conn <- data.frame(name=names(comp$membership), component=comp$membership, stringsAsFactors=FALSE)
conn$fs <- as.character(unlist(lapply(strsplit(conn$name, " "), function(l){l[1]})))
conn$date <- as.Date(unlist(lapply(strsplit(conn$name, " "), function(l){l[2]})))
conn$season <- unlist(lapply(conn$date, getSeason))
conn$cas <- apply(conn, 1, function(l){
    return(menin.cas[which(menin.cas$nomfus==as.character(l[["fs"]])),
                     l[['date']]])})
conn$pop <- apply(conn, 1, function(l){
    return(menin.pop[which(menin.pop$nomfus==as.character(l[["fs"]])),
                     l[['date']]])})
conn <- merge(conn, partition.areas.df, by.x=c("fs", "season"), by.y=c("nomfus", "season"), all.x=TRUE)
conn <- merge(conn, partition.centroids.df, by.x=c("fs", "season"), by.y=c("nomfus", "season"), all.x=TRUE)

save(comp, compp, conn, file="../data/foyers/160526_foyers.RData") ## comp, compp
