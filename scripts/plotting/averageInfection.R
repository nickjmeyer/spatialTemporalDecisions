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

data(county.fips)

avgInfObs = t(read.table("../../data/wns/obsData.txt",header=FALSE))/2
avgInfObs = as.data.frame(avgInfObs)
names(avgInfObs) = paste("Observed.Year.",as.character(2006+0:7),sep="")

avgInfSpatial = t(read.table("../../data/wns/2016-11-03-16-24-26/sampStats_mean_avgInf_gravity2_spatial_.txt",header=FALSE))
avgInfSpatial = as.data.frame(avgInfSpatial)
names(avgInfSpatial) = paste("Spatial.Year.", as.character(2006+0:7),sep="")

avgInfEdge = t(read.table("../../data/wns/2016-11-03-16-24-26/sampStats_mean_avgInf_edgeToEdge2_edge_.txt",header=FALSE))
avgInfEdge = as.data.frame(avgInfEdge)
names(avgInfEdge) = paste("Edge.Year.", as.character(2006+0:7),sep="")

avgInfOosSpatial = t(read.table("../../data/wns/2016-11-03-16-24-26/sampStats_mean_Oos_avgInf_gravity2_spatial_.txt",header=FALSE))
avgInfOosSpatial = as.data.frame(avgInfOosSpatial)
names(avgInfOosSpatial) = paste("Spatial.Oos.Year.", as.character(2012+0:1),sep="")

avgInfOosEdge = t(read.table("../../data/wns/2016-11-03-16-24-26/sampStats_mean_Oos_avgInf_edgeToEdge2_edge_.txt",header=FALSE))
avgInfOosEdge = as.data.frame(avgInfOosEdge)
names(avgInfOosEdge) = paste("Edge.Oos.Year.", as.character(2012+0:1),sep="")

avgInf = cbind(avgInfObs,avgInfSpatial,avgInfEdge,avgInfOosSpatial,avgInfOosEdge)

fips = c(as.matrix(read.table("../../data/wns/fips.txt",header=FALSE)))

d = data.frame(cbind(avgInf,fips))


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

dWNS = dWNS[order(dWNS$group,dWNS$order),]
dNA = dNA[order(dWNS$group,dWNS$order),]

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

## average infection observed
pAvgInf = ggplot()
pAvgInf = pAvgInf + geom_polygon(data=usaData, aes(x=long,y=lat,group=group,
                         order=order),color="black",fill="gray80")
pAvgInf = pAvgInf + geom_polygon(data=dWNS, aes(x=long,y=lat,group=group,
                         fill=as.factor(Observed.Year.2013),order=order),color="black")
pAvgInf = pAvgInf + settings
pAvgInf = pAvgInf +  scale_fill_manual(values=c("dodgerblue4","green"),
                                          na.value="transparent",
                                          name = "Probability of Infection\n")
pAvgInf = pAvgInf + ggtitle("Observed Average Infection")
print(pAvgInf)
ggsave("../../data/plotting/average_observed_col.pdf",plot=pAvgInf,width=width,height=height)


## average infection spatial
pAvgInf = ggplot()
pAvgInf = pAvgInf + geom_polygon(data=usaData, aes(x=long,y=lat,group=group,
                         order=order),color="black",fill="gray80")
pAvgInf = pAvgInf + geom_polygon(data=dWNS, aes(x=long,y=lat,group=group,
                         fill=log(Spatial.Year.2013+0.001),order=order),color="black")
pAvgInf = pAvgInf + settings
pAvgInf = pAvgInf +  scale_fill_gradientn(colours=c("dodgerblue4","firebrick","green"),
                                          values=rescale(
                                                   c(min(log(dWNS$Spatial.Year.2013+0.001)),
                                                     quantile(log(dWNS$Spatial.Year.2013+0.001),
                                                              probs=0.75),
                                                     max(log(dWNS$Spatial.Year.2013)))),
                                          na.value="transparent",
                                          name = "Log-Probability of Infection\n")
pAvgInf = pAvgInf + ggtitle("Spatial Average Infection")
print(pAvgInf)
ggsave("../../data/plotting/average_infection_spatial_col.pdf",plot=pAvgInf,width=width,height=height)

## average infection edge
pAvgInf = ggplot()
pAvgInf = pAvgInf + geom_polygon(data=usaData, aes(x=long,y=lat,group=group,
                         order=order),color="black",fill="gray80")
pAvgInf = pAvgInf + geom_polygon(data=dWNS, aes(x=long,y=lat,group=group,
                         fill=log(Edge.Year.2013+0.001),order=order),color="black")
pAvgInf = pAvgInf + settings
pAvgInf = pAvgInf +  scale_fill_gradientn(colours=c("dodgerblue4","firebrick","green"),
                                          values=rescale(
                                                   c(min(log(dWNS$Edge.Year.2013+0.001)),
                                                     quantile(log(dWNS$Edge.Year.2013+0.001),
                                                              probs=0.75),
                                                     max(log(dWNS$Edge.Year.2013)))),
                                          na.value="transparent",
                                          name = "Log-Probability of Infection\n")
pAvgInf = pAvgInf + ggtitle("Edge Average Infection")
print(pAvgInf)
ggsave("../../data/plotting/average_infection_edge_col.pdf",plot=pAvgInf,width=width,height=height)

## average infection spatial out-of-sample
pAvgInf = ggplot()
pAvgInf = pAvgInf + geom_polygon(data=usaData, aes(x=long,y=lat,group=group,
                         order=order),color="black",fill="gray80")
pAvgInf = pAvgInf + geom_polygon(data=dWNS, aes(x=long,y=lat,group=group,
                         fill=log(Spatial.Oos.Year.2013+0.001),order=order),color="black")
pAvgInf = pAvgInf + settings
pAvgInf = pAvgInf +  scale_fill_gradientn(colours=c("dodgerblue4","firebrick","green"),
                                          values=rescale(
                                                   c(min(log(dWNS$Spatial.Oos.Year.2013+0.001)),
                                                     quantile(log(dWNS$Spatial.Oos.Year.2013+0.001),
                                                              probs=0.5),
                                                     max(log(dWNS$Spatial.Oos.Year.2013)))),
                                          na.value="transparent",
                                          name = "Log-Probability of Infection\n")
pAvgInf = pAvgInf + ggtitle("Spatial Out-of-Sample Average Infection")
plot(pAvgInf)
ggsave("../../data/plotting/average_infection_spatial_oos_col.pdf",plot=pAvgInf,width=width,height=height)

## average infection edge out-of-sample
pAvgInf = ggplot()
pAvgInf = pAvgInf + geom_polygon(data=usaData, aes(x=long,y=lat,group=group,
                         order=order),color="black",fill="gray80")
pAvgInf = pAvgInf + geom_polygon(data=dWNS, aes(x=long,y=lat,group=group,
                         fill=log(Edge.Oos.Year.2013+0.001),order=order),color="black")
pAvgInf = pAvgInf + settings
pAvgInf = pAvgInf +  scale_fill_gradientn(colours=c("dodgerblue4","firebrick","green"),
                                          values=rescale(
                                                   c(min(log(dWNS$Edge.Oos.Year.2013+0.001)),
                                                     quantile(log(dWNS$Edge.Oos.Year.2013+0.001),
                                                              probs=0.75),
                                                     max(log(dWNS$Edge.Oos.Year.2013)))),
                                          na.value="transparent",
                                          name = "Log-Probability of Infection\n")
pAvgInf = pAvgInf + ggtitle("Edge Out-of-Sample Average Infection")
print(pAvgInf)
ggsave("../../data/plotting/average_infection_edge_oos_col.pdf",plot=pAvgInf,width=width,height=height)
