rm(list=ls(all=TRUE))

source("makeToyNets.R")

netSize = as.numeric(commandArgs(TRUE)[1])
generateAndSaveNets(netSize)
