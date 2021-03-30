## Generate Figure 13 of the paper: locations of HC
## By MW, GPLv3+, Mar. 2017
## The panels are the following:
##+ A map with the names of the districts and a dot for each HC.

##
## ==== Loading data
##
library(rgdal)

## HC locations
locali <- read.csv("../data/geolocalisation/150413_localisations_fs.csv") # This file contains 457 lines (HC) with their date of creation and last reported case.

## Shapefiles
basename <- "../data/carto"
pays <- readOGR(dsn = paste(basename, "BF", sep="/"), "country")
regions <- readOGR(dsn = paste(basename, "BF", sep="/"), "ADMIN1")
districts <- readOGR(dsn = paste(basename, "WEB", sep="/"), "BUF") # Contient les communes et les régions
districts$NAME <- tolower(districts$ADM4)
districts$NAME <- gsub("é", "e", districts$NAME)
districts$NAME <- gsub("ô", "o", districts$NAME)
districts$NAME <- gsub(" ", "_", districts$NAME)

##
## ==== Fig. 13 (location of clusters)
##
d.y <- districts$NAME %in% unique(locali$district)
dd.y <- subset(districts, d.y)
dd.n <- subset(districts, !d.y)

## Annotations
arro <- list()
arro[["OUAHIGOUYA"]] <- c(-3.18325921976955, -2.67189799004276, 14.6207301350553, 14.181795246023) # x1, x2, y1, y2
arro[["TITAO"]] <- c(-1.81962927383144, -2.14504096547576, 14.8142579666016, 14.2120666176804)
arro[["GOURCY"]] <- c(-4.08201532050149, -2.57892322100152, 14.0153027019072,   13.0920258663564)
arro[["BOULSA"]] <- c(-0.0686044568881803, -0.517982507254148, 15.0440069706061, 13.8790815294489)
arro[["SEGUENEGA"]] <- c(-1.15331009570259, -1.78863768415102, 15.1169716831756, 13.3796038971018)
arro[["YAKO"]] <- c(-4.43841860182622, -2.73388116940358, 13.3796038971018, 12.8952619505833)
arro[["DANDE"]] <- c(-5.08924198511486, -4.60887234506848, 12.562276862351, 12.047663544176)
arro[["ORODARA"]] <- c(-5.41465367675919, -5.15122516447569, 12.1687490308056, 11.9265780575464)
arro[["HOUNDÉ"]] <- c(-1.60268814606855, -2.88883911780564,10.6551804479355, 11.2908792527409)
arro[["LENA"]] <- c(-1.91260404287267, -3.8650741927386, 10.216245558903, 11.2606078810835)
arro[["KARANGASSO VIGUE"]] <- c(-2.16053676031596, -3.81858680821799, 9.82271772735688, 10.8216729920512)
arro[["DÔ"]] <- c(-5.58510742000145, -4.73283870379013, 9.9589388998152, 10.7914016203938)
arro[["DAFRA"]] <- c( -5.24419993351692, -4.20598167922313, 9.65622518324116, 10.9730298503382)
arro[["DÉDOUGOU"]] <- c(-5.15252, -3.718928, 13.26406, 12.30374)

arro.n <- names(arro)
arro <- as.data.frame(do.call("rbind", arro))
names(arro) <- c("x1", "x2", "y1", "y2")
arro$district <- arro.n
off <- c(rep(.1,8), rep(-.1,5), .1)

## Actual plot
pdf("../papers/woringer2/figures/map_hc.pdf", height=10, width=12)

plot(pays, lwd=3)
plot(dd.n, add=T, col="gray")
plot(dd.y, add=T)
plot(pays, lwd=3, add=T)
points(latitude~longitude, data=locali, cex=0.5)
arrows(arro$x1, arro$y1, arro$x2, arro$y2, length=0.1)
text(arro$x1, arro$y1+off, arro$district)

dev.off()
