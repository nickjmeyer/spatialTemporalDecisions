rm(list=ls(all=TRUE))

source("makeToyNets.R")

netSize = as.numeric(commandArgs(TRUE)[1])
nets = generateNets(netSize,FALSE)
for(i in nets)
  save(i)
