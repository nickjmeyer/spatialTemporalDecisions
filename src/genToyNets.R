rm(list=ls(all=TRUE))

source("makeToyNets.R")

netSizes = as.numeric(commandArgs(TRUE))
res = generateAndSaveNets(netSizes)
