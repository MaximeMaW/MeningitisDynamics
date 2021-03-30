## Generate Figure 7 of the paper: analysis of clusturing
## By MW, GPLv3+, Oct. 2016
## The panels are the following:
## - i.   Distribution of cluster sizes
## - ii.  Area covered vs. duration
## - iii. 10 examples of dynamics

##
## ==== Load stuff
##
load(file="../data/foyers/160517_foyers.RData") ## comp, compp
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg

##
## ==== Various functions
##
getSeason <- function(da){
    ## Return the season of a date
    y <- as.numeric(substr(da, 1, 4))
    if (da < as.Date(paste(y, 07, 01, sep="-"))) {
        return(y)
    } else {
        return(y+1)
    }
}


##
## ==== Panel 1
##
sumextent <- as.data.frame(table(compp$sumextent))
sumextent$Var1 <- as.numeric(as.character(sumextent$Var1))
sumextent$Freq <- as.numeric(sumextent$Freq)

## Prepare the plot
sum.toplot <- data.frame(Var1=1:max(sumextent$Var1))
sum.toplot <- merge(sum.toplot, sumextent, all.x=T)

## Do the plot (we split the stuff in several subplots so that we will
##+be able to work with it in Inkscape
pdf("../papers/woringer2/figures/cluster_analysis/distribution_cluster_sizes.pdf", width=10, height=6)

barplot(sum.toplot$Freq, names.arg=sum.toplot$Var1,
     xlab="Number of fs",
     ylab="Number of clusters",
     main="Distribution of the number of fs per cluster",
     col="#0093AF",
     xlim=c(0,22),
     ylim=c(0,5)) ## Number of epidemic fs

barplot(sum.toplot$Freq, names.arg=sum.toplot$Var1,
     xlab="Number of fs",
     ylab="Number of clusters",
     main="Distribution of the number of fs per cluster",
     col="#0093AF",     
     xlim=c(0,22)+80,
     ylim=c(0,5)) ## Number of epidemic fs

barplot(sum.toplot$Freq, names.arg=sum.toplot$Var1,
     xlab="Number of fs",
     ylab="Number of clusters",
     main="Distribution of the number of fs per cluster",
     col="#0093AF",     
     xlim=c(0,22),
     ylim=c(0,5)+95) ## Number of epidemic fs

dev.off()

##
## ==== Panel 2
## The idea is to plot the duration of clusters (in weeks) vs the
##+duration of the event (in weeks).

## Compute the areas of the shapefiles
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

## Preformat the connected components object
conn <- data.frame(name=names(comp$membership), component=comp$membership, stringsAsFactors=FALSE)
conn$fs <- unlist(lapply(strsplit(conn$name, " "), function(l){l[1]}))
conn$date <- as.Date(unlist(lapply(strsplit(conn$name, " "), function(l){l[2]})))
conn$season <- unlist(lapply(conn$date, getSeason))

## Merge with bobby
conn.areas <- merge(conn, partition.areas.df, by.x=c("fs", "season"), by.y=c("nomfus", "season"))

## Get the mean HC size (2012)
hc2012.mean <- mean(partition.areas[["2012"]]$area)
hc2012.median <- median(partition.areas[["2012"]]$area)
hc2012.max <- max(partition.areas[["2012"]]$area)
hc2012.min <- min(partition.areas[["2012"]]$area)

## Compute aggregated statistics
compp$maxarea <- apply(compp, 1, function(l) {
    co <- subset(conn.areas, component==as.numeric(l[['component']]))
    co.un <- co[!duplicated(co$fs),]
    return(sum(co.un$area))
})

## Make a beautiful plot (?)

pdf("../papers/woringer2/figures/cluster_analysis/duration_vs_extension.pdf", width=10, height=6)
compp$durationw <- compp$duration/7
cols.cuts <- c(1,2,3,4,5,6,7,8,9,10,15,20,80)
cols.maps <- cut(subset(compp$sumextent, compp$sumextent!=1), cols.cuts)
cols.hm <- heat.colors(13)
cols <- cols.hm[as.numeric(cols.maps)]

plot(maxarea ~ durationw, data=subset(compp, sumextent==1),
     pch=17,
     cex=1.5,
     xlim=c(0,12),
     ylim=c(0,5000)
)
points(maxarea ~ durationw, data=subset(compp, sumextent!=1),
       col=cols,
       pch=16,
       cex=2)
legend("topleft", legend=levels(cols.maps), col=cols.hm, pch=16)
legend("topright", legend=1, pch=17)
lines(c(0,6), c(hc2012.median, hc2012.median), col='black', lty=2, lwd=3)
dev.off()

##
## ==== Panel 3
## Show some of the clusters
## 

## SEE CLUSTANALIS.R
## ALL THE PANEL IS GENERATED THERE

