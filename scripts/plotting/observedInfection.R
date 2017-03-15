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
library(viridis)

data(county.fips)

obsData = t(read.table("../../data/wns/obsData.txt",header=FALSE))/2
obsData = as.data.frame(obsData)
names(obsData) = paste("Observed.Year.",as.character(2006+0:7),sep="")

yearInf = as.vector(ncol(obsData) - rowSums(as.matrix(obsData)) + 2006)
yearInf = ifelse(yearInf < 2014, yearInf, Inf)
yearInf = factor(yearInf, levels = c(2006 + 0:7, Inf),
                 labels = c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "Uninfected"))

fips = c(as.matrix(read.table("../../data/wns/fips.txt",header=FALSE)))

d = data.frame(yearInf = yearInf, fips = fips)


usaData = map_data("usa")
for(i in names(d))
    usaData[[i]] = as.numeric(rep(NA,nrow(usaData)))

fipsRegion = data.frame(matrix(unlist(strsplit(as.character(county.fips$polyname),split=",")),
    byrow=TRUE,ncol=2))
names(fipsRegion) = c("region","subregion")
fipsRegion$fips = county.fips$fips
fipsRegion$polyname = county.fips$polyname

d = merge(fipsRegion,d,by="fips")

countyData = map_data("county")

d = merge(d,countyData,by=c("region","subregion"),all=TRUE,sort=FALSE)
d = data.frame(d,polyname=paste(d$region,d$subregion,sep=","))


dWNS = d[which(!is.na(d$fips)),]
dNA = d[which(is.na(d$fips)),]
dWNSInf = d[which(!is.na(d$fips) & d$yearInf != "Uninfected"),]
dWNSNot = d[which(!is.na(d$fips) & d$yearInf == "Uninfected"),]

dWNS = dWNS[order(dWNS$group,dWNS$order),]
dNA = dNA[order(dWNS$group,dWNS$order),]
dWNSInf = dWNSInf[order(dWNSInf$group,dWNSInf$order),]
dWNSNot = dWNSNot[order(dWNSNot$group,dWNSNot$order),]



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

settings = theme(axis.line=element_blank(), axis.text.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks=element_blank(),
    axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.background=element_blank(), panel.border=element_blank(),
    panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    legend.position=c(.89,.2),
    ## legend.background=element_rect(colour="black"),
    legend.text=element_text(size=20),
    legend.key.size=unit(.5,"in"),
    legend.title=element_text(size=15))

width=20.16
height=11.13
alpha=.35
alphaCol=.35
cex=10

pObsInf = ggplot()
pObsInf = pObsInf + geom_polygon(data=dWNSNot,
                                 aes(x=long, y=lat, group=group),
                                 color="gray50",fill=NA)
pObsInf = pObsInf + geom_polygon(data=usaData,
                                 aes(x=long, y=lat, group=group),
                                 color="black",fill=NA)
pObsInf = pObsInf + geom_polygon(data=dWNSInf,
                                 aes(x=long, y=lat, group=group, fill=(as.numeric(yearInf)+2005)),
                                 color="gray50")
pObsInf = pObsInf + scale_fill_viridis(name="Year Infected")
pObsInf = pObsInf + settings
print(pObsInf)

ggsave("../../data/plotting/observed_infection_col.pdf",plot=pObsInf,width=width,height=height)
