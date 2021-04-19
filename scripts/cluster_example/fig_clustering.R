## Generate Figure 4 of the paper: clusturing.
## By MW, GPLv3+, Oct. 2016
## The panels are the following:
## - i.   Spatio temporal K-Ripley

##
## ==== Load libraries
##
library(rgdal)
library(rgeos)
library(maptools)
library(zoo)
library(fields)
library(spatstat)
library(stpp)  
source("utilitaires.R")

##
## ==== Load data
##
basename <- "../data/carto" ## Load spatial data
districts <- readOGR(dsn = paste(basename, "WEB", sep="/"), "BUF", p4s="+init=epsg:4326") # Contient les communes et les régions
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg

amh <- read.csv("../data/epi/150413_bdd_FRFE.csv") ## Load and prepare epidemiological data
amh$cas[amh$district=="titao" & amh$date=="2011-03-27"] <- 0 ## Correction d'un bug pesant
menin <- subset(amh, maladie=="menin")
fus <- read.csv("../data/geolocalisation/150413_fusion_fs_2500m.csv") ## Carte des fusions
menin.cas <- fusionner.fs(menin, fus, "cas", verbose=FALSE)
menin.pop <- fusionner.fs(menin, fus, "population", verbose=FALSE)
menin.inc <- menin.cas
idx <- 1:(ncol(menin.inc)-1)
menin.inc[,idx] <- menin.cas[,idx]/menin.pop[,idx]

geo <- aggregate(fus$latitude, by=list(fus$nomfus), FUN=mean) ## Aggregate spatial data
lo <- aggregate(fus$longitude, by=list(fus$nomfus), FUN=mean)
names(geo) <- c('nomfus', 'latitude')
geo$longitude <-  lo$x


districts.rep <- spTransform(districts, CRS("+init=epsg:32630")) ## Reprojections and co.
districts.rep <- gBuffer(districts.rep, width=0, byid=T) # Debug shapefile
geo2 <- SpatialPointsDataFrame(geo[,c("longitude","latitude")], geo, proj4string=CRS("+init=epsg:4326")) # 
geo.proj <- spTransform(geo2, CRS("+init=epsg:32630"))

##
## ==== Functions
##
## Détermine si les semaines sont épidémiques
is.epi.week.seuiladj <- function(cas, pop, qq) {
    ## Découpage en saisons
    y.ref <- as.Date(names(cas)[1:(length(cas)-1)])
    saisons <- list()
    for (a in 2004:2015) {
        idx1 <- (y.ref < as.Date(paste(a,"07", "01", sep="-"))) & (y.ref > as.Date(paste(a-1,"06", "30", sep="-")))
        saisons[[as.character(a)]] <- cas[,c(idx1,TRUE)]
        saisons[[as.character(a)]]$saison <- a
        saisons[[as.character(a)]]$population <- pop[,which.max(idx1)]
    }
    ## Seuillage
    seu <- lapply(saisons, function(d) {e <- d
                                        e$popseuil <- qbinom(0.95, round(e$population), qq)
                                        idxneg <- -1*c(ncol(e)-1, ncol(e)-2, ncol(e)-3)
                                        res <- t(apply(e[,idxneg], 1, function(w)
                                            {return(as.numeric(w[1:(length(w)-1)]) >= w[length(w)])}))
                                        e[,1:(ncol(e)-4)] <- res
                                        return(e)
                                    })
    ## Reconstruction
    nomfus <-seu[[1]]$nomfus
    seu2 <- lapply(seu, function(e) {return(e[,1:(ncol(e)-4)])})
    cn <- as.character(unlist(lapply(seu2, names)))
    out <- do.call("cbind", seu2)
    names(out) <- cn
    out$nomfus <- nomfus

    return(out)
}

##
## ==== Threshold epidemiological data
##
qq <- 5/1e5
qq40 <- 40/1e5
seu <- is.epi.week.seuiladj(menin.cas,menin.pop,qq)
seu40 <- is.epi.week.seuiladj(menin.cas,menin.pop,qq40)
seu30 <- is.epi.week.seuiladj(menin.cas,menin.pop,30/1e5)


## -- Interesting representation (but irrelevant here)
## seupl <- as.matrix(seu[,1:(ncol(seu)-1)])
## rownames(seupl) <- seu$nomfus
## cn <- substr(colnames(seupl), 1,7)
## par(mfrow=c(1,2))
## image(seupl, axes=F, main=paste("Semaines épidémiques, seuil : ", qq), col=rev(heat.colors(128)))
## mtext(text=cn[seq(1,ncol(seupl),by=20)],
##       side=4, line=0.3, at=seq(0,1,20/ncol(seupl)), las=1, cex=0.7)
## image.plot(seupl, main=paste("Semaines épidémiques, seuil : ", qq), axes=F)

##
## ===== Compute spatio-temporal K-Ripley
##

## Prepare epidemiological data
stp <- melt(seu30) # flatten the array
stp <- na.omit(stp) # remove NA
stp <- stp[stp$value==TRUE,1:2] # keep epidemic weeks
names(stp) <- c('nomfus', 'date') # rename
fusu <- as.data.frame(coordinates(geo.proj))
fusu$nomfus <- geo.proj$nomfus
stpm <- merge(stp, fusu[,c('nomfus', 'latitude', 'longitude')], by='nomfus') # merge

xyt <- stpm[,c('longitude', 'latitude', 'date')]
names(xyt) <- c('x','y','t')
xytpp <- as.3dpoints(xyt) # How to incorporate bounding box information?

##plot(xytpp)
## Compute bounding polygon -- avoid singular point pattern
FMD <- as.3dpoints(xyt$x/1000, xyt$y/1000, xyt$t)
offx <- rnorm(nrow(xyt), 1/1000)
offy <- rnorm(nrow(xyt), 1/1000)
xyt$x <- xyt$x+offx ## Avoid overlapping points
xyt$y <- xyt$y+offy
FMD <- as.3dpoints(xyt$x/1000, xyt$y/1000, as.Date(xyt$t))
ttt <- as.numeric(as.Date(FMD[,3]))

## Select the districts for each group
##n1 <- c("Boulsa")
n2 <- c("Ouahigouya", "Titao", "Yako", "Gourcy", "Seguenega")
N2 <- c("ouahigouya", "titao", "yako", "gourcy", "seguenega")
n2.name <- "North"
n3 <- c("Orodara", "Dédougou", "Karangasso Vigué", "Léna", "Houndé", "Dô", "Dafra", "Dandé")
N3 <- c("orodara", "dedougou", "karangasso vigue", "lena", "hounde", "do", "dafra", "dande")
n3.name <- "West"

##s1 <- subset(districts.rep, ADM4 %in% n1)
s2 <- subset(districts.rep, ADM4 %in% n2)
s3 <- subset(districts.rep, ADM4 %in% n3)
##ss1 <- unionSpatialPolygons(s1, s1$ADM4 %in% n1)
ss2 <- unionSpatialPolygons(s2, s2$ADM4 %in% n2)
ss3 <- unionSpatialPolygons(s3, s3$ADM4 %in% n3)
##sp1 <- ss1@polygons[[1]]@Polygons[[1]]@coords
sp2 <- ss2@polygons[[1]]@Polygons[[1]]@coords
sp3 <- ss3@polygons[[1]]@Polygons[[1]]@coords
##ppp1 <- ppp(xyt$x, xyt$y,window=owin(poly=sp1[rev(1:nrow(sp1)),]), marks=ttt)
ppp2 <- ppp(xyt$x, xyt$y,window=owin(poly=sp2[rev(1:nrow(sp2)),]), marks=ttt)
ppp3 <- ppp(xyt$x, xyt$y,window=owin(poly=sp3[rev(1:nrow(sp3)),]), marks=ttt)
##pp1 <- cbind(coords(ppp1)/1000, marks(ppp1))
pp2 <- cbind(coords(ppp2)/1000, marks(ppp2))
pp3 <- cbind(coords(ppp3)/1000, marks(ppp3))

## ==== Compute the K-Ripley & Resample
## Generate new stuff

## Compute the K-Ripley
u <- seq(0, 250, by = 1) 
v <- seq(0, 6*7, by = 7) # times
nrep <- 20

##stik1 <- STIKhat(xyt = pp1,
##                 dist = u, times = v, infectious = TRUE,
##                 s.region=sp1[rev(1:nrow(sp1)),]/1000)
stik2 <- STIKhat(xyt = pp2,
                 dist = u, times = v, infectious = TRUE,
                 s.region=sp2[rev(1:nrow(sp2)),]/1000)
stik3 <- STIKhat(xyt = pp3,
                 dist = u, times = v, infectious = TRUE,
                 s.region=sp3[rev(1:nrow(sp3)),]/1000)

## Commented because this is only Boulsa (very small number of points)
## pprandom1 <- list()
## for (i in 1:nrep) {
##     print(paste("-- iteration", i, "/", nrep))
##     pp1p <- coords(rpoispp(.1, win=owin(poly=sp1[rev(1:nrow(sp1)),]/1000))[1:nrow(pp1)])
##     pp1p$t <- pp1[,3]
##     pprandom1[[i]] <- STIKhat(xyt = pp1p,
##                               dist = u, times = v, infectious = TRUE,
##                               s.region=sp1[rev(1:nrow(sp1)),]/1000)
## }

resampleEpidemics <- function(ppnp, pp.res) {
    ## Function resample a dataframe, reasmpling is performed on a timepoint basis
    ## Resampling is performed uniformly among the HC locations.
    pp <- as.data.frame(table(ppnp[,3])) ## compute the number of
    res <- apply(pp, 1, function(l){
        res <- pp.res[sample(1:nrow(pp.res), as.numeric(l["Freq"])),]
        res$t <- as.numeric(as.character(l["Var1"]))
        return(res)
    })
    ress <- do.call("rbind", res)
    names(ress) <- c("x","y","t")
    return(ress)
}

getSeason <- function(da){
    ## Return the season of a date
    y <- as.numeric(substr(da, 1, 4))
    if (da < as.Date(paste(y, 07, 01, sep="-"))) {
        return(y)
    } else {
        return(y+1)
    }
}
extractAdjacencyInit <-  function() {
    for (a in as.character(2004:2014)) {
        saisons.adj[[a]] <- gTouches(partition.merg[[a]], byid=TRUE)
    }
}
geo.proj.df <- as.data.frame(coordinates(geo.proj))/1000
geo.proj.df$nomfus <- geo.proj$nomfus
names(geo.proj.df) <- c("x", "y", "nomfus")
saisons.adj <- list()
extractAdjacencyInit() ## Init saisons.adj

reverseLookup <- function(xx, yy){
    subs <- subset(geo.proj.df,
                   abs(geo.proj.df$x-as.numeric(xx))<.1 & abs(geo.proj.df$y-as.numeric(yy))<.1)
    if (nrow(subs)==1) {
        return(subs)
    } else {
        error("CPT")
    }
}

extractAdjacency <- function(xx,yy,d) {
    ## Function returns a dataframe of adjacent locations
    xx <- as.numeric(xx)
    yy <- as.numeric(yy)
    d <- as.Date(d)
    
    ## Reverse lookup the location
    subs <- subset(geo.proj.df,
                   abs(geo.proj.df$x-as.numeric(xx))<.1 & abs(geo.proj.df$y-as.numeric(yy))<.1)
    if (nrow(subs)==1) {
         ## Find the corresponding adjacency HC
        y <- substr(d, 1,4)
        fs <- as.character(subs$nomfus)
        co <- which(colnames(saisons.adj[[y]])==fs) ## column to extract
        adj.names <- c(colnames(saisons.adj[[y]])[as.logical(saisons.adj[[y]][co,])], fs)

        ## Forward lookup
        adj <- subset(geo.proj.df, nomfus %in% adj.names)
        return(adj)
    } else {
        error("CPT")
    }   
}

resampleEpidemicsTC <- function(ppnp, pp.res, tc=0, sc=0) {
    ## Function resample a dataframe, reasmpling is performed on a timepoint basis
    ## Resampling is performed with a specified time correlation among the HC
    ## locations.
    pp <- as.data.frame(table(ppnp[,3])) ## compute the number of HC per timepoint
    res <- list()
    for (i in 1:nrow(pp)) {
        l <- pp[i,]
        t.cur <- as.numeric(as.character(pp$Var1[i]))
        fr <- as.numeric(l["Freq"])
        if (i>1 && t.cur-t.old==7) {
            inh <- res[[i-1]][runif(nrow(res[[i-1]]))<tc,1:2]
            if (nrow(inh)>= fr) { # We have enough or too much
                res[[i]] <- inh[sample(1:nrow(inh), fr),] ## Save the line
            } else { ## We need more (do not add duplicates)
                fr.m <- as.integer((fr-nrow(inh))*(1-sc)) ## How many points to put at random
                fr.n <- fr-fr.m-nrow(inh) ## How many points to correlate
                dups <- apply(pp.res, 1, function(ll) {
                    if ((ll["x"] %in% inh$x) && (ll["y"] %in% inh$y)) {
                        return(FALSE)
                    } else {
                        return(TRUE)
                    }
                })
                pp.new <- pp.res[dups,]
                inh.new <- pp.new[sample(1:nrow(pp.new), fr.m),]/1000.

                ## Iterative approach to add the correlated points
                tmp <- rbind(inh, inh.new)
                if (nrow(tmp)==0) { ## Seed the algorithm
                    tmp <- pp.new[sample(1:nrow(pp.new), 1),]
                }
                tmp$t <- as.Date(t.cur)
                tmp.loc <- apply(tmp, 1, function(l){
                    reverseLookup(as.numeric(l["x"]),
                                  as.numeric(l["y"]))})
                tmp.loc <- do.call("rbind", tmp.loc)
                if (fr.n>0) {
                    for (k in 1:fr.n) {
                        ## Extract unique adjacency
                        ok <- FALSE
                        while (!ok) {
                            adj <- do.call("rbind", apply(tmp, 1, function(l){
                                extractAdjacency(as.numeric(l["x"]),
                                                 as.numeric(l["y"]),
                                                 t.cur)
                            })) ## Extract all the adjacent stuff
                            adj <- adj[!duplicated(adj$nomfus) & !duplicated(adj$nomfus, fromLast=T),] ## Remove duplicates
                            dups <- as.logical(
                                apply(adj, 1, function(l){
                                    l["nomfus"] %in% tmp.loc$nomfus}))
                            adj <- adj[!dups,] ## Remove more duplicates
                            if (nrow(adj)!=0) {
                                ok <- TRUE
                            } else {
                                print("BAD, reseeding adjacency")
                                tmp2 <- pp.new[sample(1:nrow(pp.new), 1),]
                                tmp2$t <- as.Date(t.cur)
                                tmp2.loc <- reverseLookup(
                                    as.numeric(tmp2["x"]),
                                    as.numeric(tmp2["y"]))
                                tmp <- rbind(tmp, tmp2)
                                tmp.loc <- rbind(tmp.loc, tmp2.loc)
                            }
                        }
                        names(adj) <- c("x", "y", "nomfus")
                        adj[,1:2] <- adj[,1:2]*1000
                        adj$t <- as.Date(t.cur)
                        sel <- adj[sample(1:nrow(adj),1),]
                        sel2 <- sel
                        sel2[,1:2] <- sel2[,1:2]/1000
                        names(sel2) <- c("x", "y", "nomfus", "t")
                        tmp <- rbind(tmp, sel[c("x", "y", "t")])
                        tmp.loc <- rbind(tmp.loc, c(sel2[c("x", "y", "nomfus")]))
                    }
                    tmp.loc <- tmp.loc[sample(1:nrow(tmp.loc), fr.n),]
                    tmp <- tmp[sample(1:nrow(tmp), fr.n),]
                }
                inh.sc <- tmp.loc[,c("x", "y")]
                inh.sc$t <- as.numeric(tmp$t)
                res[[i]] <- inh.sc ## Save the line
            }
        } else {
            res[[i]] <- pp.res[sample(1:nrow(pp.res), fr),]/1000
            names(res[[i]]) <- c("x", "y")
        }
        res[[i]]$t <- t.cur
        t.old <- t.cur
    }
    ress <- do.call("rbind", res)
    ress[,1:2] <- ress[,1:2] + matrix(rnorm(2*nrow(ress), sd=1/1000), ncol=2)
    names(ress) <- c("x","y","t")
    return(ress)
}

## Collect the estimates
plotresample <- function(pprandom, nn, stik, bar=NULL, xlim=c(0,150)) {
    utr <- stik$dist
    Khats <- list()
    Khatsflat <- list()
    for (i in 1:nrep) {
        Khats[[i]] <- pprandom[[i]]$Khat
        Khatsflat[[i]] <- matrix(pprandom[[i]]$Khat, nrow=1)
    }
    Khatstab <- do.call('rbind', Khatsflat)
    Khatsquantflat <- apply(Khatstab, 2, function(l){quantile(l, probs=.95)})
    Khatsquantflatlow <- apply(Khatstab, 2, function(l){quantile(l, probs=.05)})
    Khatsquant <- matrix(Khatsquantflat, nrow=nrow(Khats[[1]]))
    Khatsquantlow <- matrix(Khatsquantflatlow, nrow=nrow(Khats[[1]]))
    
    tshift <- c(1,2,3)  ## This is weeks!!

    plot(range(utr), c(0,max(stik$Khat)),
         pch='', main=nn, xlab='distance (km)', ylab='cumulated #pairs',
         xlim=xlim)
    for (i in 1:length(tshift)) {
        j <- tshift[i]
        polygon(c(utr, rev(utr)),
                c(Khatsquant[,j], rev(Khatsquantlow[,j])),
                col = rgb(0,0,0,20,maxColorValue=255), border = NA)
        lines(utr, stik$Khat[,j], lty=i) # -stik$Khat[1,j])
        lines(utr, (Khatsquant[,j]), lty=i, col='red')
        lines(utr, (Khatsquantlow[,j]), lty=i, col='red')
    }
    lines(range(utr), c(1,1))
    if (!is.null(bar)) {
        lines(c(bar, bar), c(0, max(stik$Khat)*1.5), col="blue", lwd=2)
    }

    legend('topleft', lty=rep(1:length(tshift),2),
           col=rep(c('black', 'red'), each=length(tshift)),
           legend=c(paste('Delta t (week):',tshift), paste('Ref. Delta t (week):',tshift)))
}    

plotresamplediff <- function(pprandom, nn, stik, bar=NULL, xlim=c(0,150)) {
    utr <- stik$dist
    Khats <- list()
    Khatsflat <- list()
    for (i in 1:nrep) {
        Khats[[i]] <- pprandom[[i]]$Khat
        Khatsflat[[i]] <- matrix(pprandom[[i]]$Khat, nrow=1)
    }
    Khatstab <- do.call('rbind', Khatsflat)
    Khatsquantflat <- apply(Khatstab, 2, function(l){quantile(l, probs=.95)})
    Khatsquantflatlow <- apply(Khatstab, 2, function(l){quantile(l, probs=.05)})
    Khatsquant <- matrix(Khatsquantflat, nrow=nrow(Khats[[1]]))
    Khatsquantlow <- matrix(Khatsquantflatlow, nrow=nrow(Khats[[1]]))
    
    tshift <- c(1,2,3)  ## This is weeks!!

    plot(range(utr), c(0,max(diff(stik$Khat))),
         pch='', main=nn, xlab='distance (km)', ylab='cumulated #pairs',
         xlim=xlim)
    for (i in 1:length(tshift)) {
        j <- tshift[i]
        polygon(c(utr[1:(length(utr)-1)], rev(utr[1:(length(utr)-1)])),
                c(diff(Khatsquant[,j]), rev(diff(Khatsquantlow[,j]))),
                col = rgb(0,0,0,20,maxColorValue=255), border = NA)
        lines(utr[1:(length(utr)-1)], diff(stik$Khat[,j]), lty=i) # -stik$Khat[1,j])
        lines(utr[1:(length(utr)-1)], diff(Khatsquant[,j]), lty=i, col='red')
        lines(utr[1:(length(utr)-1)], diff(Khatsquantlow[,j]), lty=i, col='red')
    }
    lines(range(utr), c(1,1))
    if (!is.null(bar)) {
        lines(c(bar, bar), c(0, max(stik$Khat)*1.5), col="blue", lwd=2)
    }

    legend('topleft', lty=rep(1:length(tshift),2),
           col=rep(c('black', 'red'), each=length(tshift)),
           legend=c(paste('Delta t (week):',tshift), paste('Ref. Delta t (week):',tshift)))
}  

resampleEpidemicsBourrin <- function(pp2, pp2.res, s=10) {
    res <- apply(as.data.frame(table(pp2[,3])), 1, function(l) {
        ##sa <- as.numeric(pp2.res[sample(1:nrow(pp2.res), 1),])/1000
        sa <- colMeans(pp2.res)/1000
        re <- as.data.frame(matrix(rep(sa,
                                       as.numeric(l["Freq"])),
                                   ncol=2, byrow=TRUE))
        re <- re+matrix(rnorm(2*nrow(re), sd=s), ncol=2)
        re$t <- as.numeric(l["Var1"])
        return(re)
    })
    return(do.call("rbind", res))
}

resampleEpidemicsSC <-function(pp2, pp2.res, sp2, s=10) {
    # Resample epidemics for the pure SC model (no time correlation)
    res <- apply(as.data.frame(table(pp2[,3])), 1, function(l) {
        ##sa <- as.numeric(pp2.res[sample(1:nrow(pp2.res), 1),])/1000
        p = runifpoint(n=1, win=owin(poly=sp2[rev(1:nrow(sp2)),]))
        sa <- c(p$x,p$y)/1000
        re <- as.data.frame(matrix(rep(sa,
                                       as.numeric(l["Freq"])),
                                   ncol=2, byrow=TRUE))
        re <- re+matrix(rnorm(2*nrow(re), sd=s), ncol=2)
        re$t <- as.numeric(l["Var1"])
        return(re)
    })
    return(do.call("rbind", res))
}    


resampleEpidemicsMulti <- function(pp2, pp2.res, s=c(10,20), p=0.5) {
    ## Function that generates a mix of two types of clusters of characteristic
    ## size described in variable `s`, and a proportion described by `p`.
    n1 <- rbinom(1, nrow(pp2), p)
    n2 <- nrow(pp2)-n1
    
    pp2s1 <- pp2[1:n1,]
    pp2s2 <- pp2[(n1+1):nrow(pp2),]

    if (n2 != 0) {
        return(rbind(resampleEpidemicsBourrin(pp2s1, pp2.res, s[1]),
                     resampleEpidemicsBourrin(pp2s2, pp2.res, s[2])))
    } else {
        resampleEpidemicsBourrin(pp2s1, pp2.res, s[1])
    }
}


saveResample <- FALSE

if (saveResample) {
    ## Get ready for the resampling
    geo.proj$district <- unlist(lapply(strsplit(as.character(geo.proj$nomfus), "__"), function(l){l[[1]]}))
    
    pprandom2 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemics(pp2, pp2.res) ## Custom function, see above
        pp2p[,c("x","y")] <- pp2p[,c("x","y")]/1000+matrix(rnorm(2*nrow(pp2p), sd=1/1000), ncol=2)

        pprandom2[[i]] <- STIKhat(xyt = pp2p,
                                  dist = u, times = v, infectious = TRUE,
                                  s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }

    pprandom2.11 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsBourrin(pp2, pp2.res, s=11) ## Custom function, see above
        pprandom2.11[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom2.10 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsBourrin(pp2, pp2.res, s=10) ## Custom function, see above
        pprandom2.10[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom2.9 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsBourrin(pp2, pp2.res, s=9) ## Custom function, see above
        pprandom2.9[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom2.8.11.p0.5 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsMulti(pp2, pp2.res, s=c(8,11), p=0.5)
        pprandom2.8.11.p0.5[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom2.11.13.p0.5 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsMulti(pp2, pp2.res, s=c(11,13), p=0.5)
        pprandom2.11.13.p0.5[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }    
    pprandom2.10.13.p0.15 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsMulti(pp2, pp2.res, s=c(10,13), p=0.15)
        pprandom2.10.13.p0.15[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom2.10.13.p0.85 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsMulti(pp2, pp2.res, s=c(10,13), p=0.85)
        pprandom2.10.13.p0.85[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom2.8.13.p0.85 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        pp2p <- resampleEpidemicsMulti(pp2, pp2.res, s=c(8,13), p=0.85)
        pprandom2.8.13.p0.85[[i]] <- STIKhat(xyt = pp2p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    pprandom3.20 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        pp3p <- resampleEpidemicsBourrin(pp3, pp3.res, s=20) ## Custom function, see above
        pprandom3.20[[i]] <- STIKhat(xyt = pp3p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }
    pprandom3.25 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        pp3p <- resampleEpidemicsBourrin(pp3, pp3.res, s=25) ## Custom function, see above
        pprandom3.25[[i]] <- STIKhat(xyt = pp3p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }    
    pprandom3.30 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        pp3p <- resampleEpidemicsBourrin(pp3, pp3.res, s=30)
        pprandom3.30[[i]] <- STIKhat(xyt = pp3p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }    
    pprandom3.35 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        pp3p <- resampleEpidemicsBourrin(pp3, pp3.res, s=35) ## Custom function, see above
        pprandom3.35[[i]] <- STIKhat(xyt = pp3p,
                                     dist = u, times = v, infectious = TRUE,
                                     s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }
    
    pprandom3 <- list()
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        pp3p <- resampleEpidemics(pp3, pp3.res) ## Custom function, see above
        pp3p[,c("x","y")] <- pp3p[,c("x","y")]/1000+matrix(rnorm(2*nrow(pp3p), sd=.1), ncol=2)
        
        pprandom3[[i]] <- STIKhat(xyt = pp3p,
                                  dist = u, times = v, infectious = TRUE,
                                  s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }
    #save(n2, n3, stik2, stik3, nrep,
    #     pprandom2, pprandom3,
    #     pprandom2.8, pprandom2.9, pprandom2.10, pprandom2.11, pprandom2.14, 
    #     pprandom3.20, pprandom3.25, pprandom3.35, pprandom3.35,
    #     pprandom2.8.13.p0.85, pprandom2.10.13.p0.85, pprandom2.10.13.p0.15,
    #     file="../data/clusters/170302_stKripley_resample.RData")
    save(n2, n3, stik2, stik3, nrep, 
         pprandom2, pprandom2.10, pprandom2.10.13.p0.15, pprandom2.10.13.p0.85, pprandom2.11, 
         pprandom2.11.13.p0.5, pprandom2.8.11.p0.5, pprandom2.8.13.p0.85, 
         pprandom2.9, pprandom3, pprandom3.20, pprandom3.25, pprandom3.30, pprandom3.35, 
         file="../data/clusters/210414_stKripley_resample.RData")
} else {
    ## n2, n3, stik2, stik3, pprandom2, pprandom3,
    load("../data/clusters/210414_stKripley_resample.RData")
}

## Addition 2021
if (saveResample) {
    ## Get ready for the resampling
    geo.proj$district <- unlist(lapply(strsplit(as.character(geo.proj$nomfus), "__"), function(l){l[[1]]}))

    pprandom2.tc0 <- list()
    print("pprandom2.tc0")
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        names(pp2.res)<-c('x','y')
        pp2p <- resampleEpidemicsTC(pp2, pp2.res, tc=0) ## Custom function, see above
        pp2p[,c("x","y")] <- pp2p[,c("x","y")]+matrix(rnorm(2*nrow(pp2p), sd=1/1000), ncol=2)
        
        pprandom2.tc0[[i]] <- STIKhat(xyt = pp2p,
                                      dist = u, times = v, infectious = TRUE,
                                      s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    
    pprandom2.tc1 <- list()
    print("pprandom2.tc1")
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        names(pp2.res)<-c('x','y')
        pp2p <- resampleEpidemicsTC(pp2, pp2.res, tc=1) ## Custom function, see above
        pp2p[,c("x","y")] <- pp2p[,c("x","y")]+matrix(rnorm(2*nrow(pp2p), sd=1/1000), ncol=2)
        
        pprandom2.tc1[[i]] <- STIKhat(xyt = pp2p,
                                      dist = u, times = v, infectious = TRUE,
                                      s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
    
    pprandom3.tc0 <- list()
    print("pprandom3.tc0")
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        names(pp3.res)<-c('x','y')
        pp3p <- resampleEpidemicsTC(pp3, pp3.res, tc=0) ## Custom function, see above
        pp3p[,c("x","y")] <- pp3p[,c("x","y")]+matrix(rnorm(2*nrow(pp3p), sd=.1), ncol=2)
        
        pprandom3.tc0[[i]] <- STIKhat(xyt = pp3p,
                                      dist = u, times = v, infectious = TRUE,
                                      s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }

    pprandom3.sc30 <- list()
    print("pprandom3.sc30")
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        names(pp3.res)<-c('x','y')
        pp3p <- resampleEpidemicsSC(pp3, pp3.res, sp3, s=30) ## Custom function, see above
        names(pp3p)<-c('x','y','t')
        pp3p[,c("x","y")] <- pp3p[,c("x","y")]+matrix(rnorm(2*nrow(pp3p), sd=.1), ncol=2)
        
        pprandom3.sc30[[i]] <- STIKhat(xyt = pp3p,
                                      dist = u, times = v, infectious = TRUE,
                                      s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }    

    pprandom2.sc10 <- list()
    print("pprandom2.sc10")
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp2.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N2)))
        names(pp2.res)<-c('x','y')
        pp2p <- resampleEpidemicsSC(pp2, pp2.res, sp2, s=10) ## Custom function, see above
        names(pp2p)<-c('x','y','t')
        pp2p[,c("x","y")] <- pp2p[,c("x","y")]+matrix(rnorm(2*nrow(pp2p), sd=.1), ncol=2)
        
        pprandom2.sc10[[i]] <- STIKhat(xyt = pp2p,
                                       dist = u, times = v, infectious = TRUE,
                                       s.region=sp2[rev(1:nrow(sp2)),]/1000)
    }
            
    pprandom3.tc1 <- list()
    print("pprandom3.tc1")
    for (i in 1:nrep) {
        print(paste("-- iteration", i, "/", nrep))
        pp3.res <- as.data.frame(coordinates(subset(geo.proj, district %in% N3)))
        names(pp3.res)<-c('x','y')
        pp3p <- resampleEpidemicsTC(pp3, pp3.res, tc=1) ## Custom function, see above
        pp3p[,c("x","y")] <- pp3p[,c("x","y")]+matrix(rnorm(2*nrow(pp3p), sd=.1), ncol=2)
        
        pprandom3.tc1[[i]] <- STIKhat(xyt = pp3p,
                                      dist = u, times = v, infectious = TRUE,
                                      s.region=sp3[rev(1:nrow(sp3)),]/1000)
    }
    save(n2, n3, stik2, stik3, nrep, 
         pprandom2.tc0, pprandom2.tc1,
         pprandom3.tc0, pprandom3.tc0, 
         pprandom2.sc10, pprandom3.sc30,
         file="../data/clusters/210414_stKripley_resampleTC.RData")
} else {
    ## n2, n3, stik2, stik3, pprandom2, pprandom3,
    load("../data/clusters/210414_stKripley_resampleTC.RData")
}

##
## ==== Plot stuff
##
plotresample(pprandom2.8.13.p0.85, n2, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.10.13.p0.15, n2, stik2, bar=70, xlim=c(0,100))

##
## ==== Figure 6
##
pdf("../papers/woringer2/figures/cluster_example/kripley-simple-raw.pdf", height=6, width=13)
par(mfrow=c(1,2))
plotresample(pprandom3.30, n3.name, stik3, bar=NA, xlim=c(0,150)) ## West
plotresample(pprandom2.10, n2.name, stik2, bar=NA, xlim=c(0,100)) ## North
dev.off()

pdf("../papers/woringer2/figures/cluster_example/kripley.pdf", height=8, width=12)
par(mfrow=c(1,2))
plotresample(pprandom2, n2, stik2, bar=70, xlim=c(0,100)) 
plotresample(pprandom3, n3, stik3, bar=140)
dev.off()

## Do the resampling
pdf("../figures/cluster_example/Kripley_zone2.pdf", height=8, width=24)
par(mfrow=c(2,3))
plotresample(pprandom2,    n2, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.8.11.p0.5,    n2, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.9, n2, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.11, n2, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.8.13.p0.85, n2, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.10.13.p0.85, n2, stik2, bar=70, xlim=c(0,100))
#plotresample(pprandom2.14.11.p25, n2, stik2, bar=70, xlim=c(0,100))
dev.off()

pdf("../figures/cluster_example/Kripley_zone3.pdf", height=8, width=16)
par(mfrow=c(1,2))
plotresample(pprandom3, n3, stik3, bar=140, xlim=c(0,150))
plotresample(pprandom3.30, n3, stik3, bar=140, xlim=c(0,150))
dev.off()

pdf("../figures/cluster_example/Kripley_vs_null_raw.pdf", height=4, width=8)
par(mfrow=c(1,2))
plotresample(pprandom2, n2.name, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom3, n3.name, stik3, bar=140, xlim=c(0,150))
dev.off()

pdf("../figures/cluster_example/Kripley_STC_model.pdf", height=8, width=24)
par(mfrow=c(1,2))
plotresample(pprandom2.11, n2.name, stik2, bar=70, xlim=c(0,100))
#plotresample(pprandom2.14, n2.name, stik2, bar=70, xlim=c(0,100))
plotresample(pprandom3.35, n3.name, stik3, bar=140, xlim=c(0,150))
dev.off()


## Same with empirical H curve
##par(mfrow=c(1,2))
##plotresamplediff(pprandom2, n2, stik2, bar=70, xlim=c(0,100))
##plotresamplediff(pprandom3, n3, stik3, bar=140)

##
## ==== First try with a simple model for temporal aggregation
##


pdf("../figures/cluster_example/Kripley_STC_model.pdf", height=12, width=28)
par(mfrow=c(2,4))
plotresample(pprandom2.tc0, "CSR - North", stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.tc1, "TC - North", stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.sc10, "SC - North", stik2, bar=70, xlim=c(0,100))
plotresample(pprandom2.10, "STC, 10km clusters - North", stik2, bar=70, xlim=c(0,100))

plotresample(pprandom3.tc0, "CSR - West", stik3, bar=140, xlim=c(0,150))
plotresample(pprandom3.tc1, "TC - West", stik3, bar=140, xlim=c(0,150))
plotresample(pprandom3.sc30, "SC - West", stik3, bar=140, xlim=c(0,150))
plotresample(pprandom3.30, "STC, 30km clusters - West", stik3, bar=140, xlim=c(0,150))
dev.off()


#stik2.tc1 <- STIKhat(xyt = resampleEpidemicsTC(pp2, pp2.res, tc=1),
#                     dist = u, times = v, infectious = TRUE,
#                     s.region=sp2[rev(1:nrow(sp2)),]/1000)



#stik2.tc0 <- STIKhat(xyt = resampleEpidemicsTC(pp2, pp2.res, tc=0),
#                  dist = u, times = v, infectious = TRUE,
#                  s.region=sp2[rev(1:nrow(sp2)),]/1000
                  )
## stik2.sc1 <- STIKhat(xyt = resampleEpidemicsTC(pp2, pp2.res, tc=0, sc=.1),
##                  dist = u, times = v, infectious = TRUE,
##                  s.region=sp2[rev(1:nrow(sp2)),]/1000)
## stik2.bourr1 <- STIKhat(xyt = resampleEpidemicsBourrin(pp2, pp2.res, s=14),
##                  dist = u, times = v, infectious = TRUE,
##                  s.region=sp2[rev(1:nrow(sp2)),]/1000)
## stik2.bourr2 <- STIKhat(xyt = resampleEpidemicsBourrin(pp2, pp2.res, s=11),
##                  dist = u, times = v, infectious = TRUE,
##                  s.region=sp2[rev(1:nrow(sp2)),]/1000)
## stik3.bourr1 <- STIKhat(xyt = resampleEpidemicsBourrin(pp3, pp3.res, s=20),
##                  dist = u, times = v, infectious = TRUE,
##                  s.region=sp3[rev(1:nrow(sp3)),]/1000)


## Quiche/quick compare
## par(mfrow=c(1,2))
## plotresample(pprandom2, n2, stik2, bar=70, xlim=c(0,100))
## lines(stik2.bourr1$dist, stik2.bourr1$Khat[,1], col='orange')
## lines(stik2.bourr1$dist, stik2.bourr1$Khat[,2], col='orange')
## lines(stik2.bourr1$dist, stik2.bourr1$Khat[,3], col='orange')
## lines(stik2.bourr2$dist, stik2.bourr2$Khat[,1], col='yellow')
## lines(stik2.bourr2$dist, stik2.bourr2$Khat[,2], col='yellow')
## lines(stik2.bourr2$dist, stik2.bourr2$Khat[,3], col='yellow')

## plotresample(pprandom3, n3, stik3, bar=140, xlim=c(0,150))
## lines(stik3.bourr1$dist, stik3.bourr1$Khat[,1], col='orange')
## lines(stik3.bourr1$dist, stik3.bourr1$Khat[,2], col='orange')
## lines(stik3.bourr1$dist, stik3.bourr1$Khat[,3], col='orange')

## plotresample(pprandom2, n2, stik2.bourr, bar=70, xlim=c(0,100))
## plotresample(pprandom2, n2, stik2.tc1, bar=70, xlim=c(0,100))
## plotresample(pprandom2, n2, stik2.tc0, bar=70, xlim=c(0,100))
## plotresample(pprandom2, n2, stik2.sc1, bar=70, xlim=c(0,100))

## Conclusion: temporal aggregation do not explain the enrichment pattern we observe.
## Conclusion: a simple model of (very strong) spatial aggregation explains the stuff
