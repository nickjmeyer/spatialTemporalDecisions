rm(list=ls(all=TRUE))


generateNets <- function(n,display=TRUE){
  nets = c("alleyNet","bowTieNet","gridNet","randNet",
           "scaleFreeNet")
  nets = foreach(net = nets)%dopar%{
    source("makeToyNets.R")

    argGen = get(paste(net,"Args",sep=""))
    netGen = get(paste("gen",
                       paste(toupper(substring(net,1,1)),
                             substring(net,2),sep=""),sep=""))

    netRes = do.call(netGen,argGen(n))
    cat(paste(net,"is done","\n"))
    return(netRes)
  }

  if(display){
    for(net in nets){
      plotNet(net)
      readline(prompt="Press [enter] to continue.")
    }
    dev.off()
  }

  return(nets)
}



netSize = as.numeric(commandArgs(TRUE)[1])
nets = generateNets(netSize,FALSE)
for(i in nets)
  saveNet(i)
