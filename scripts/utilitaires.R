## Utilitaires pour bosser avec la base de données
## By MW, GPLv3+
## Apr. 2015

## Chargement des librairies
library(reshape)
library(zoo)
library(RColorBrewer)
library(classInt)


## Découpage en saisons
to.seasons <- function(menin.cas, miny=2004, maxy=2015) {
    y.ref <- as.Date(names(menin.cas)[1:(length(menin.cas)-1)])
    saisons <- list()
    for (a in miny:maxy) {
        idx1 <- (y.ref < as.Date(paste(a,"07", "01", sep="-"))) & (y.ref > as.Date(paste(a-1,"06", "30", sep="-")))
        saisons[[as.character(a)]] <- menin.cas[,c(idx1,TRUE)]
        saisons[[as.character(a)]]$annee <- a
    }
    return(saisons)
}

## Détermine si une semaine est épidémique en regardant par rapport à un seuil ajusté
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

## Usage 
#qq <- 5/1e5
#seu <- is.epi.week.seuiladj(cas,pop,qq)
#seupl <- as.matrix(seu[,1:(ncol(seu)-1)])
#rownames(seupl) <- seu$nomfus


## Calculer si une semaine est épidémique ou non. Cette fonction est un peu pourrie parce
## qu'elle ne prend pas finement en compte le fait que la population évolue au cours du temps
is.epi.week.seuiladj.old <- function(cas, pop, qq) {
    out <- pop
    nc <- ncol(cas)-1
    for (i in 1:nrow(cas)) {
        ccas <- cas[i,1:nc]
        ppop <- pop[i,1:nc]
        sseuil <- qbinom(0.95, as.numeric(round(ppop)), qq)
        out[i, 1:(length(out)-1)] <- ccas > sseuil
    }
    return(out)
}


## Récupérer les coordonnées des séries temporelles
## Moyenner quand il y en a plusieurs
## Entrées :
## - le nom des formations sanitaires fusionnées (nom standard)
## - le tableau "fus" qui permet de revenir en arrière sur la fusion
## Sorties :
## - Un dataframe c(latitude, longitude)
## - où les coordonnées sont la moyenne des coordonnées des points fusionnés
coordonnees.fus <- function(noms, fus) {
    ## Fonction qui retourne un dataframe c("noms", "latitude", "longitude")
    ## Avec les coordonnées correspondantes au blob fusionné "noms"
    r <- lapply(noms, function(n) {
        colMeans(subset(fus, nomfus==n, c("latitude", "longitude")))
    })
    r <- as.data.frame(do.call("rbind", r))
    cbind(noms, r)
}


## Appliquer la fusion des formations sanitaires
## Entrées :
##   - le tableau des fusions
##   - un export de la base de données
## Sortie
##   - un tableau de la base de données : nom, date, cas
sum2 <- function(l) { ## Fusionner bien comme il faut
    if (all(is.na(l))) {
        return(NA)
    } else {
        return(sum(l, na.rm=TRUE))
    }
}


fusionner.fs <- function(dat, fus, col="cas", verbose=TRUE) {
    ## Créer la matrice au format zoo (série temporelle)
    fus$nomfus <- as.character(fus$nomfus)
    ts.dat <- cast(dat[,c("district", "fsjoli", "date", col)], district+fsjoli~date)
    rownames(ts.dat) <- paste(ts.dat$district, ts.dat$fsjoli, sep="__")
    ts.dat$fus <- paste(ts.dat$district, ts.dat$fsjoli, sep=".")

    ## Fusionner les formations sanitaires
    ts.out <- ts.dat[0,3:(ncol(ts.dat)-1)]
    ts.out$nomfus <- character(0)
    i <- 1
    for (l in unique(fus$nomfus)) { ## Boucler sur les noms fusionnés
        ss <- subset(fus, nomfus==l)
        sub <- subset(ts.dat, fus %in% ss$fus) ## Un subset des séries temporelles
        sub2 <- subset(sub, TRUE, 3:(ncol(sub)-1)) ## Prendre les bonnes colonnes
        ll <- c(apply(sub2, 2, sum2), l) ## Sommer
        ts.out[i,] <- ll ## Injecter dans le tableau
        if (verbose) {
            print(paste(i, ":", l, nrow(sub)))
        }
        i <- i+1
        ## tester pour voir si on a le bon nombre de lignes à chaque fois, sinon, quelque chose a merdé
    }
    ts.out[,1:(ncol(ts.out)-1)] <- apply(ts.out[,1:(ncol(ts.out)-1)], 2, as.numeric)
    return(ts.out)
}

## Afficher une carte pour une date donnée
## Penser à pouvoir régler et afficher une légende comme il faut
carte <- function(ti, df.m, shp,
                  maladie="maladie non renseignée",
                  max.scale, breaks=c(0,1,2,3,5,7,9,11,15,20,30,50,60,80,100),
                  legend.intervals=TRUE, legend.seuil=FALSE) {
    
    if (!(as.character(ti) %in% names(df))) {
        return("La date sélectionnée n'existe pas")
    }
    
    ts.ti <- subset(df.m, date==as.Date(ti), c("nomfus", "cas", "col")) #df[,c("nomfus",ti)]
        
    names(ts.ti) <- c("NAME", "cas", "col")
    y <- substr(as.character(ti), 1, 4)
    shp.ti <- shp[[y]]
    shp.ti$cas <- merge(shp.ti, ts.ti, by="NAME", all.x=TRUE)$cas
    shp.ti$col <- merge(shp.ti, ts.ti, by="NAME", all.x=TRUE)$col
    #shp.ti$col <- rev(heat.colors(max(shp.ti$cas, na.rm=TRUE)+1))[shp.ti$cas+1]
    
    plot(shp.ti, col=shp.ti$col, main=paste(maladie, ti))
    plot(districts, lwd=2, add=TRUE)
    if (legend.intervals) {
        legend('topleft', legend=c(names(attr(colcode, 'table')),'no data'), fill=c(attr(colcode, 'palette'),'white'), cex=.7)
    } else if (legend.seuil==TRUE) {
        legend("topleft", legend=c("épidémique", "non épidémique", "non renseigné"), col=c("red", "white", "grey"), pch=15)
    }
    #image.plot(legend.only=TRUE, zlim=range(shp.ti$cas, na.rm=T), col=attr(colcode, "palette"), horizontal=T)
}
