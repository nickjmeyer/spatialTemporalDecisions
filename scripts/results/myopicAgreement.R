rm(list=ls(all=TRUE))

## dataDir = "../../data/wns/2017-01-26-08-40-46/"
dataDir = "../../data/wns/2017-01-26-09-29-40/"


myopicData = as.matrix(read.table(paste(dataDir,"myopic_trt_.txt",sep="")))
rankData = as.matrix(read.table(paste(dataDir,"rank_trt_.txt",sep="")))

myopicInf = myopicData >= 2
rankInf = rankData >= 2

myopicTrt = myopicData %% 2 == 1
rankTrt = rankData %% 2 == 1

mean(myopicTrt == rankTrt)

mean((myopicTrt == rankTrt)[which(rankInf)])
mean((myopicTrt == rankTrt)[which(!rankInf)])
