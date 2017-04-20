rm(list=ls(all=TRUE))

library(ggplot2)
library(maps)
library(maptools)
library(sp)
library(scales)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(gganimate)

data(county.fips)

dataDir = "../.././data/wns/./2017-03-20-11-15-53/"

files = dir(dataDir)

repNum = 0

rankData = read.table(
             paste(dataDir, "rank_M1Sp_history_", repNum, "_.txt",
                   sep = ""))

fips = c(as.matrix(read.table("../../data/wns/fips.txt",header=FALSE)))

names(rankData) = fips

rankData$Year = 2006 + 1:nrow(rankData) - 1

rankData = melt(rankData, id.vars = c("Year"),
     variable.name = "fips",
     value.name = "Status")

rankData$Status = factor(rankData$Status, levels = 0:3,
                         labels = c("Uninfected", "Uninfected w/ Trt",
                                    "Infected", "Infected w/ Trt"))

usaData = map_data("usa")

fipsRegion = data.frame(matrix(unlist(strsplit(as.character(county.fips$polyname),split=",")),
    byrow=TRUE,ncol=2))
names(fipsRegion) = c("region","subregion")
fipsRegion$fips = county.fips$fips
fipsRegion$polyname = county.fips$polyname

dWNS = merge(fipsRegion,rankData,by="fips", all.y = TRUE, sort=FALSE)

countyData = map_data("county")

dWNS = merge(dWNS,countyData,by=c("region","subregion"),all.x=TRUE,sort=FALSE)
dWNS = data.frame(dWNS,polyname=paste(dWNS$region,dWNS$subregion,sep=","))

dWNS = dWNS[order(dWNS$fips, dWNS$group, dWNS$order), ]


getLabelPoint = function(county){Polygon(county[c("long","lat")])@labpt}
countyCentroids = dWNS[which(!is.na(dWNS$fips) & !is.na(dWNS$lat)
    & !is.na(dWNS$long)),]
countyCentroids = countyCentroids[order(countyCentroids$group,countyCentroids$order),]
countyCentroids = by(countyCentroids,countyCentroids$polyname,getLabelPoint)
countyCentroids = do.call("rbind.data.frame", countyCentroids)
names(countyCentroids) = c("long","lat")

## p = qplot(long,lat,data=d,geom="polygon",fill=value,group=group,order=order)
## p = p + geom_polygon(data=dVal, aes(x=long,y=lat,group=group,order=order),color="black")
## p = p + geom_polygon(data=dNA, aes(x=long,y=lat,group=group,order=order))
## p = p + geom_polygon(data=usaData, aes(x=long,y=lat,group=group,order=order),color="black")
## p = p + theme_bw()

width=20.16
height=11.13
alpha=.35
alphaCol=.35
cex=10

settings = theme(axis.line=element_blank(), axis.text.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks=element_blank(),
    axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.background=element_blank(), panel.border=element_blank(),
    panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    legend.position=c(.89,.2),
    ## legend.background=element_rect(colour="black"),
    legend.text=element_text(size=16),
    legend.key.size=unit(.5,"in"),
    legend.title=element_text(size=17),
    plot.title = element_text(hjust = 0.5, size = 18))

p = ggplot()
p = p + geom_polygon(data=usaData, aes(x=long,y=lat,group=group),
                     color="black",fill="gray80")
p = p + geom_polygon(data=dWNS, aes(x=long,y=lat,group=group,
                         fill=Status,frame=Year),color="black")
p = p + settings
p = p + ggtitle("Year ")


animation::ani.options(ani.width = width * 50, ani.height = height * 50)
gganimate(p, "../../data/animations/animate_spread.gif")
