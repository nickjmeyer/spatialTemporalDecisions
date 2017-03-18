rm(list=ls(all=TRUE))
library(doMC)
library(Matrix)
library(MASS)
library(ggplot2)
library(igraph)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(truncnorm)
registerDoMC(detectCores())


sourceCpp("./main/getCov.cpp")
sourceCpp("./main/mergeClusters.cpp")

genRandNet<-function(n,numNeigh=3){
  set.seed(0)

  nodes=matrix(runif(2*n),ncol=2)
  d=as.matrix(dist(nodes))
  neigh=diag(n)
  for(i in 1:n){
    ind=order(d[i,])[2:(numNeigh+1)]
    neigh[i,ind]=1
    neigh[ind,i]=1
  }



  clusters = getClusters(neigh,n)
  nClusters = length(unique(clusters))

  while(nClusters > 1){
    neigh = mergeOneCluster(neigh,clusters,nodes[,1],nodes[,2],n)
    neigh = matrix(neigh,ncol=n)

    clusters = getClusters(neigh,n)
    nClusters = length(unique(clusters))
  }


  cov=getCov(n,nodes)

  net=list(n=n,neigh=neigh,nodes=nodes,
           fips=1:n,caves=cov$caves,Xcov=cov$Xcov)
  net$gDist=getDist(net)

  meanCaves=mean(net$caves)

  net$start = getStart(net$n)

  net$name = "rand"

  net$hasNetwork = TRUE;

  return(net)
}



genAlleyNetBreak<-function(nodes){

  net=genAlleyNet(nodes)

  n=net$n

  if(n==100)
    netDir="starredAlleyBreakSmall"
  else if(n==500)
    netDir="starredAlleyBreakMedium"
  else if(n==1000)
    netDir="starredAlleyBreakLarge"
  else
    stop("invalid n")

  files=dir(netDir)
  if("wnsSubGraph.txt" %in% files){
    subGraph=read.table(paste(netDir,"/wnsSubGraph.txt",sep=""),header=FALSE)
    between=read.table(paste(netDir,"/wnsBetweenness.txt",sep=""),header=FALSE)

    subGraph=max(subGraph)-subGraph
    between=max(between)-between

    subGraph=((subGraph - min(subGraph))/(max(subGraph)-min(subGraph)))^4
    between=((between - min(between))/(max(between)-min(between)))^4

    subGraph=ifelse(subGraph<.5,0,subGraph)
    between=ifelse(between<.5,0,between)

    measure=colMeans(rbind(c(subGraph),c(between)))

    caves=20*measure + 1
    caves=round(caves)
  }
  else{
    caves=rep(1,n)
    for(i in 1:n)
      caves[i]=round(10*sum((abs(nodes[i,]-c(.5,.25)))^c(1.0,.25)))+1
    cat("No sub graph measure file.  Run again to generate the correct caves.\n")
  }

  net$caves=caves
  net$Xcov[,1]=scale(caves)

  net$start = getStart(net$n)

  net$hasNetwork = TRUE;

  return(net)
}



genRingNet<-function(n1,n2){
  n=n1+n2
  theta=2*pi/(n1+1)
  radius = 1/theta ## arc length = 1
  rotMat=getRotMat(theta)
  nodes=c(0,radius)
  curNode=nodes
  for(i in 2:n1){
    curNode=c(rotMat%*%curNode)
    nodes=rbind(nodes,curNode)
  }

  theta=theta/n2
  rotMat=getRotMat(theta)
  for(i in 1:n2){
    curNode=c(rotMat%*%curNode)
    nodes=rbind(nodes,curNode)
  }

  neigh=diag(n)
  for(i in 1:n){
    if(i>1){
      neigh[i,i-1]=1
      neigh[i-1,i]=1
    }
    if(i<n){
      neigh[i,i+1]=1
      neigh[i+1,i]=1
    }
  }
  neigh[1,n]=1
  neigh[n,1]=1

  cov=getCov(n,nodes)

  net=list(n=n,neigh=neigh,nodes=nodes,
           fips=1:n,caves=cov$caves,Xcov=cov$Xcov)


  net$gDist = matrix(0,nrow=n,ncol=n)
  net$gDist[which(neigh==1)] = 1
  for(i in 0:(n2-1)){
    net$gDist[n-i,n-i-1] = 1/n2
    net$gDist[n-i-1,n-i] = net$gDist[n-i,n-i-1]
  }
  net$gDist[which(neigh==0)] = -1
  diag(net$gDist)=0

  net$gDist=getDist(net,preAlloc=1)

  meanCaves=mean(net$caves)
  net$gDist=net$gDist

  net$start = getStart(net$n)

  net$name="ring"

  net$hasNetwork = TRUE;

  return(net)
}



genGridNet<-function(n1,n2){
  n=n1*n2
  nodes=expand.grid(seq(0,1,length.out=n1),seq(0,1,length.out=n2))
  d=as.matrix(dist(nodes))
  neigh=diag(n)
  diag=0
  tol=.00001
  for(i in 1:(n)){
    xSame=abs(nodes[,1]-nodes[i,1])<tol
    ySame=abs(nodes[,2]-nodes[i,2])<tol
    xOne=abs(nodes[,1]-nodes[i,1])<(1/(n1-1)+tol)
    yOne=abs(nodes[,2]-nodes[i,2])<(1/(n2-1)+tol)
    ind=NULL
    ind=c(ind,which(xSame==TRUE & ySame==FALSE & yOne==TRUE))
    ind=c(ind,which(ySame==TRUE & xSame==FALSE & xOne==TRUE))
    if(diag==1)
      ind=c(ind,which(xOne==TRUE & yOne==TRUE))
    neigh[i,ind]=1
    neigh[ind,i]=1
  }

  cov=getCov(n,nodes)

  net=list(n=n,neigh=neigh,nodes=as.matrix(nodes),
           fips=1:n,caves=cov$caves,Xcov=cov$Xcov)

  net$gDist=getDist(net)

  meanCaves=mean(net$caves)
  net$gDist=net$gDist

  net$start = getStart(net$n)

  net$name="grid"

  net$hasNetwork = TRUE;

  return(net)
}



genCrpNet<-function(n) {
  set.seed(0)
  theta = 2

  expNumTablesOptim<-function(parTheta) {
    return ((parTheta * (digamma(parTheta + n) - digamma(parTheta)) - 2*log(n))**2)
  }
  expNumTablesGrad<-function(parTheta) {
    inside = parTheta * (digamma(parTheta + n) - digamma(parTheta)) - 2*log(n)
    chainFirst = parTheta * (trigamma(parTheta + n) - trigamma(parTheta))
    chainSecond = digamma(parTheta + n) - digamma(parTheta)
    return (2 * inside * (chainFirst + chainSecond))
  }

  out = optim(par=theta,fn=expNumTablesOptim,gr=expNumTablesGrad,
              lower=0.1,method="L-BFGS-B")

  stopifnot(out$convergence == 0)

  theta = out$par[1]

  nodes = matrix(0,nrow=n,ncol=2)
  groupId = rep(0,n)

  groups = list()
  ngroups = 0
  for(i in 1:n) {
    newProb = theta/(theta + i - 1)
    draw = runif(1)
    if(draw <= newProb || ngroups == 0) {
      ngroups = ngroups + 1
      groups[[ngroups]] = c(i)
      groupId[i] = ngroups

      nodes[i,] = runif(2)
    } else {
      joinProbs = c()
      for(j in 1:ngroups) {
        joinProbs = c(joinProbs,length(groups[[j]])/
                                (theta + i - 1))
      }
      joinProbs = joinProbs/sum(joinProbs)
      g = sample(1:ngroups,1,prob=joinProbs)
      groups[[g]] = c(groups[[g]],i)

      meanX = nodes[groups[[g]][1],1]
      meanY = nodes[groups[[g]][1],2]

      x = rtruncnorm(1,a=0.0,b=1.0,mean=meanX,sd=0.05)
      y = rtruncnorm(1,a=0.0,b=1.0,mean=meanY,sd=0.05)
      nodes[i,] = c(x,y)
      groupId[i] = g
    }
  }

  ## qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  ## col_vector = unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
  ## cols = sample(col_vector,ngroups)
  ## plot(nodes,col=cols[groupId],pch=17,xlim=c(0,1),ylim=c(0,1))


  cov=getCov(n,nodes)

  net=list(n=n,neigh=c(),nodes=nodes,fips=1:n,caves=cov$caves,Xcov=cov$Xcov)

  net$gDist=c()

  meanCaves=mean(net$caves)

  net$start=getStart(net$n)

  net$name="crp"

  net$hasNetwork = FALSE;

  return(net)
}


genAlleyNet<-function(nodes){
### makes entire alley way

  nextMain=c(0,0)
  nextMainInd=1
  alley=makeAlleyNode(nextMain,nodes[[1]])
  for(i in 2:length(nodes)){
    prevMainInd=nextMainInd
    nextMainInd=alley$n+1
    nextMain=nextMain+c(2,0)
    nextNode=makeAlleyNode(nextMain,nodes[[i]])


    ## combine coords
    alley$nodes=rbind(alley$nodes,nextNode$nodes)

    ## combine networkMat
    alley$neigh=as.matrix(bdiag(alley$neigh,nextNode$neigh))
    alley$neigh[prevMainInd,nextMainInd]=1
    alley$neigh[nextMainInd,prevMainInd]=1

    ## combine numNodes
    alley$n=alley$n+nextNode$n
  }

  cov=getCov(alley$n,alley$nodes)

  net=alley
  net$caves=cov$caves
  net$fips=1:alley$n
  net$Xcov=cov$Xcov

  net$gDist = matrix(-1,nrow=net$n,ncol=net$n)
  net$gDist[which(net$neigh == 1)] = 1
  prongs = which(abs(net$nodes[,2]) > .0000001)
  for(i in prongs){
    for(j in 1:net$n){
      if(net$neigh[i,j] == 1)
        net$gDist[i,j] = net$gDist[j,i] = .9
    }
  }
  diag(net$gDist) = 0
  net$gDist=getDist(net,preAlloc=1)

  meanCaves=mean(net$caves)
  net$gDist=net$gDist

  net$start = getStart(net$n)

  net$name="alley"

  net$hasNetwork = TRUE;

  return(net)

}



makeAlleyNode<-function(mainNode,subAlleyLens){
### mainNode is a cartesian coordinate
### subAlleyLens is the number of alleys off of mainNode and
###		the length of each sub alley

  if(sum(subAlleyLens)==0)
    return(list(numNodes=1,nodes=mainNode,neigh=as.matrix(1)))

  numSubAlleys=length(subAlleyLens)

  neigh=diag(1+sum(subAlleyLens))

  nodes=mainNode

  if((numSubAlleys %% 2) && (numSubAlleys>1)){ ## odd && more than 1
    numBot=(numSubAlleys+1)/2
    numTop=numSubAlleys-numBot
  }
  else if(!(numSubAlleys%%2)){ ## even
    numBot=numSubAlleys/2
    numTop=numSubAlleys-numBot
  }
  else{ ## 1
    numBot=1
    numTop=0
  }

  ## get angle differences
  thetaBot=pi/(numBot+1)
  if(numTop>0)
    thetaTop=pi/(numTop+1)
  else
    thetaTop=NULL


  numNodes=1
  ## make bottom
  for(i in 1:numBot){
    numNodes=numNodes+1
    lenSubAlley=subAlleyLens[i]
    curTheta=thetaBot*i
    rotMat=getRotMat(curTheta) ## get rotation matrix
    nextNode=rotMat%*%c(-1,0)+mainNode
    nodes=rbind(nodes,c(nextNode))

    neigh[1,numNodes]=neigh[numNodes,1]=1
    if(lenSubAlley>1){
      for(j in 1:(lenSubAlley-1)){
        numNodes=numNodes+1
        nodes=rbind(nodes,c(nextNode-c(0,j)))
        neigh[numNodes,numNodes-1]=neigh[numNodes-1,numNodes]=1
      }
    }
  }

  ## make top
  if(numTop>0){
    for(i in 1:numTop){
      numNodes=numNodes+1
      lenSubAlley=subAlleyLens[i+numBot]
      curTheta=thetaTop*i
      rotMat=getRotMat(curTheta) ## get rotation matrix
      nextNode=(rotMat%*%c(1,0))+mainNode
      nodes=rbind(nodes,c(nextNode))

      neigh[1,numNodes]=neigh[numNodes,1]=1
      if(lenSubAlley>1){
        for(j in 1:(lenSubAlley-1)){
          numNodes=numNodes+1
          nodes=rbind(nodes,c(nextNode+c(0,j)))
          neigh[numNodes,numNodes-1]=neigh[numNodes-1,numNodes]=1
        }
      }
    }
  }

  return (list(n=numNodes,nodes=nodes,neigh=neigh))
}



getRotMat<-function(theta){
  return(matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),nrow=2))
}



genBowTieNet<-function(grid1N,grid2N,midN){
  grid1=genGridNet(grid1N[1],grid1N[2])
  grid2=genGridNet(grid2N[1],grid2N[2])
  mid=genRandNet(midN)

  grid1LinkNum=ceiling(grid1N[2]/3)
  grid1Link=seq(from=grid1$n-grid1N[1]*grid1LinkNum,
                by=-grid1N[1],length.out=grid1LinkNum)

  grid2LinkNum=ceiling(grid2N[2]/3)
  grid2Link=seq(from=grid2$n-grid2N[1]*(grid1LinkNum+1)+1,
                by=-grid2N[1],length.out=grid2LinkNum)

  mid1LinkNum=ceiling(midN/5)
  mid1Link=order(mid$nodes[,1],decreasing=FALSE)[1:mid1LinkNum]

  mid2LinkNum=ceiling(midN/5)
  mid2Link=order(mid$nodes[,1],decreasing=TRUE)[1:mid2LinkNum]


  ## new neighbor matrix
  neigh=as.matrix(bdiag(grid1$neigh,mid$neigh,grid2$neigh))

  link1=expand.grid(grid1Link,mid1Link)
  link2=expand.grid(mid2Link,grid2Link)

  for(i in 1:nrow(link1)){
    neigh[link1[i,1],link1[i,2]+grid1$n]=1
    neigh[link1[i,2]+grid1$n,link1[i,1]]=1
  }
  for(i in 1:nrow(link2)){
    neigh[link2[i,1]+grid1$n,link2[i,2]+grid1$n+mid$n]=1
    neigh[link2[i,2]+grid1$n+mid$n,link2[i,1]+grid1$n]=1
  }

  ## scale coordinates of middle
  midScale=2.0
  mid$nodes = mid$nodes/midScale

  ## get shifts and raises for the coordinates of the middle and second grid
  shift1=diff(range(mid$nodes[,1]))/2
  shift2=shift1

  midRaise=mean(c(grid1$nodes[grid1Link,2],grid2$nodes[grid2Link,2]))-
      mean(mid$nodes[c(mid1Link,mid2Link),2])
  midShift=shift1 - (min(mid$nodes[,1]) - max(grid1$nodes[,1]))

  grid2Shift=shift2 - (min(grid2$nodes[,1]) - max(mid$nodes[,1])) +
      midShift


  ## shift and raise coordinates
  mid$nodes[,1]=mid$nodes[,1]+midShift
  mid$nodes[,2]=mid$nodes[,2]+midRaise

  grid2$nodes[,1]=grid2$nodes[,1]+grid2Shift

  ## combine coordinates
  colnames(grid1$nodes)=colnames(mid$nodes)=colnames(grid2$nodes)=c("x","y")

  nodes=grid1$nodes
  nodes=rbind(nodes,mid$nodes)
  nodes=rbind(nodes,grid2$nodes)

  n=grid1$n+mid$n+grid2$n

  Xcov=grid1$Xcov
  Xcov=rbind(Xcov,mid$Xcov)
  Xcov=rbind(Xcov,grid2$Xcov)

  caves=grid1$caves
  caves=c(caves,mid$caves)
  caves=c(caves,grid2$caves)

  net=list(n=n,neigh=neigh,nodes=as.matrix(nodes),
           fips=1:n,caves=caves,Xcov=Xcov)

  net$gDistist = getDist(net)

  meanCaves=mean(net$caves)
  net$gDist=net$gDist

  ## net$gDistist=(getDist(net)*max(grid1N,grid2N))^3

  net$start = getStart(net$n)

  net$name="bowtie"

  net$hasNetwork = TRUE;

  return(net)
}



genScaleFreeNet<-function(n){
  set.seed(0)
  netData = barabasi.game(n)

  neigh=as.matrix(get.adjacency(netData))
  neigh=neigh+t(neigh)
  diag(neigh)=1

  set.seed(0)
  net=list(n=n,neigh=neigh,
           nodes=as.matrix(scale(layout.fruchterman.reingold(netData))),
           fips=1:n)

  cov=getCov(n,net$nodes)

  net$caves=cov$caves
  net$Xcov=cov$Xcov

  net$gDist=getDist(net)

  meanCaves=mean(net$caves)
  net$gDist=net$gDist

  net$start = getStart(net$n)

  net$name="scalefree"

  net$hasNetwork = TRUE;

  return(net)
}



genShrinkNet<-function(n1,n0=2){
  net=list()

  net$n=n1^2*n0^2

  net$nodes=foreach(i=1:n1,.combine=rbind)%do%{
    if(i>1)
      nodes1=foreach(j=1:(i-1),.combine=rbind)%do%{
        foreach(k=1:n0,.combine=rbind)%do%{
          foreach(l=1:n0,.combine=rbind)%do%{
            return(c((i-1)*n0+(k-1),(j-1)*n0+(l-1)))
          }
        }
      }
    else
      nodes1=NULL

    nodes2=foreach(j=1:i,.combine=rbind)%do%{
      foreach(k=1:n0,.combine=rbind)%do%{
        foreach(l=1:n0,.combine=rbind)%do%{
          return(c((j-1)*n0+(k-1),(i-1)*n0+(l-1)))
        }
      }
    }
    return(rbind(nodes1,nodes2))
  }

  net$gDist=as.matrix(dist(net$nodes,method="manhattan"))

  net$neigh=as.matrix(net$gDist<(1+.000001))*1

  net$fips=1:(net$n)

  cov=getCov(net$n,net$nodes)

  net$caves=cov$caves
  net$Xcov=cov$Xcov

  net$start = getStart(net$n)

  net$name="shrink"

  net$hasNetwork = TRUE;

  return(net)
}



getDist<-function(net,preAlloc=0){
  if(!preAlloc){
    dyn.load("main/getDist.so")
    out=.C("getDist",
           d=as.double(rep(0,net$n*net$n)),
           n=as.integer(net$n),
           neigh=as.integer(net$neigh),
           nodesX=as.double(net$nodes[,1]),
           nodesY=as.double(net$nodes[,2]),
           preAlloc=as.integer(preAlloc))
    dyn.unload("main/getDist.so")
  }
  else{
    dyn.load("main/getDist.so")
    out=.C("getDist",
           d=as.double(net$gDist),
           n=as.integer(net$n),
           neigh=as.integer(net$neigh),
           nodesX=as.double(net$nodes[,1]),
           nodesY=as.double(net$nodes[,2]),
           preAlloc=as.integer(preAlloc))
    dyn.unload("main/getDist.so")
  }

  return(matrix(out$d,nrow=net$n))
}


getCov<-function(n,nodes,rho=10,tau=log(10),eta=log(4),p=4,fast=TRUE){
  ## tau is set so that a distance of one results in correlation of
  ## 0.1 and eta is set so that a covariate difference of one results
  ## in correaltion of 0.25
  if(fast)
    return(getCovFast(n,nodes,rho,tau,eta,p))
  else
    return(getCovOrig(n,nodes,rho,tau,eta,p))
}


getCovOrig<-function(n,nodes,rho,tau,eta,p){
  set.seed(0);
  np=n*p

  for(i in 1:ncol(nodes)){
    if(length(unique(nodes[,i])) > 1)
      nodes[,i] = scale(nodes[,i])
  }
  nodes = as.matrix(nodes)

  drift=nodes
  drift[,1]=drift[,1]^2

  sigma=matrix(0,np,np)
  mu=c(sapply(drift%*%matrix(c(2,1),ncol=1),rep,times=p))
  for(i in 1:n){
    for(j in 1:p){
      for(k in 1:n){
        for(l in 1:p){
          ij=(i-1)*p + j
          kl=(k-1)*p + l
          sigma[ij,kl]=rho*exp(-tau*sqrt(sum((nodes[i,]-nodes[k,])^2))
                               -eta*abs(j-l))
          sigma[kl,ij]=sigma[ij,kl]
        }
      }
    }
  }

  Xcov=matrix(mvrnorm(n=1,mu=mu,Sigma=sigma),nrow=n,ncol=p,byrow=TRUE)
  caves=floor(Xcov[,1]-min(Xcov[,1]) + 1)
  for(i in 1:ncol(Xcov))
    if(length(unique(Xcov[,i]))>1)
      Xcov[,i]=scale(Xcov[,i])

  return(list(Xcov=Xcov,caves=caves))
}




getCovFast<-function(n,nodes,rho,tau,eta,p){
  set.seed(0);
  np=n*p

  for(i in 1:ncol(nodes)){
    if(length(unique(nodes[,i])) > 1)
      nodes[,i] = scale(nodes[,i])
  }
  nodes = as.matrix(nodes)


  rv = rnorm(np)

  nodesX = nodes[,1]
  nodesY = nodes[,2]

  Xcov = getCovCpp(rv,nodesX,nodesY,n,rho,tau,eta,p)

  Xcov = matrix(Xcov,ncol=p,byrow=TRUE)

  caves=floor(Xcov[,1]-min(Xcov[,1]) + 1)
  for(i in 1:ncol(Xcov))
    if(length(unique(Xcov[,i]))>1)
      Xcov[,i]=scale(Xcov[,i])

  return(list(Xcov=Xcov,caves=caves))
}



getStart<-function(n){
  set.seed(0)
  return(sample((1:n)-1,ceiling(0.01*n)))
}



saveNet<-function(net,dir=NULL){
  if(is.null(dir))
    dir=paste(net$name,net$n,sep="")
  dir=paste("../data/toy/",dir,"/",sep="")
  system(paste("mkdir -p",dir))

  ## fips
  file=paste(dir,"fips.txt",sep="")
  write.table(net$fips,file,col.names=FALSE,row.names=FALSE)

  ## network
  file=paste(dir,"network.txt",sep="")
  write.table(net$neigh,file,col.names=FALSE,row.names=FALSE)

  ## gDist
  file=paste(dir,"gDist.txt",sep="")
  write.table(net$gDist,file,col.names=FALSE,row.names=FALSE)

  ## eDist
  file=paste(dir,"eDist.txt",sep="")
  write.table(as.matrix(dist(net$nodes)),file,col.names=FALSE,row.names=FALSE)

  ## caves
  file=paste(dir,"caves.txt",sep="")
  write.table(net$caves,file,col.names=FALSE,row.names=FALSE)

  ## Xcov
  file=paste(dir,"xcov.txt",sep="")
  write.table(net$Xcov,file,col.names=FALSE,row.names=FALSE)

  ## centroids
  file=paste(dir,"centroids.txt",sep="")
  write.table(net$nodes,file,col.names=FALSE,row.names=FALSE)

  ## centroidsLong
  file=paste(dir,"centroidsLong.txt",sep="")
  write.table(net$nodes[,1],file,col.names=FALSE,row.names=FALSE)

  ## centroidsLat
  file=paste(dir,"centroidsLat.txt",sep="")
  write.table(net$nodes[,2],file,col.names=FALSE,row.names=FALSE)

  ## starting locations
  file=paste(dir,"startingLocations.txt",sep="")
  write.table(net$start,file,col.names=FALSE,row.names=FALSE)

  ## trtStart
  file=paste(dir,"trtStart.txt",sep="")
  write.table(7,file,col.names=FALSE,row.names=FALSE)

  ## period
  file=paste(dir,"period.txt",sep="")
  write.table(2,file,col.names=FALSE,row.names=FALSE)

  if(net$hasNetwork) {
    ## get betweenness connectivity
    betweenness = betweenness(graph.adjacency(net$neigh),nobigint=FALSE)
    file = paste(dir,"betweenness.txt",sep="")
    write.table(betweenness,file,col.names=FALSE,row.names=FALSE)

    ## get subGraph centrality
    subGraph = subgraph.centrality(graph.adjacency(net$neigh),diag=TRUE)
    file = paste(dir,"subGraph.txt",sep="")
    write.table(subGraph,file,col.names=FALSE,row.names=FALSE)

    ## centroidsMds
    centroidsMds = cmdscale(net$gDist,k=2)
    file=paste(dir,"centroidsMds.txt",sep="")
    write.table(centroidsMds,file,col.names=FALSE,row.names=FALSE)

    ## centroidsMdsLong
    file=paste(dir,"centroidsMdsLong.txt",sep="")
    write.table(centroidsMds[,1],file,col.names=FALSE,row.names=FALSE)

    ## centroidsMdsLat
    file=paste(dir,"centroidsMdsLat.txt",sep="")
    write.table(centroidsMds[,2],file,col.names=FALSE,row.names=FALSE)
  } else {
    ## get betweenness connectivity
    file = paste(dir,"betweenness.txt",sep="")
    write.table(c(),file,col.names=FALSE,row.names=FALSE)

    ## get subGraph centrality
    file = paste(dir,"subGraph.txt",sep="")
    write.table(c(),file,col.names=FALSE,row.names=FALSE)

    ## centroidsMds
    file=paste(dir,"centroidsMds.txt",sep="")
    write.table(c(),file,col.names=FALSE,row.names=FALSE)

    ## centroidsMdsLong
    file=paste(dir,"centroidsMdsLong.txt",sep="")
    write.table(c(),file,col.names=FALSE,row.names=FALSE)

    ## centroidsMdsLat
    file=paste(dir,"centroidsMdsLat.txt",sep="")
    write.table(c(),file,col.names=FALSE,row.names=FALSE)
  }

  ## prior treatment mean
  file = paste(dir,"priorTrtMean.txt",sep="")
  write.table(0.0,file,col.names=FALSE,row.names=FALSE)
}



plotNetNoAdj<-function(net){
  plot(net$nodes[,1],net$nodes[,2],pch=19,cex=0.5,
       xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
}



plotNet<-function(net,...){
  if (length(net$neigh) == 0) {
    plotNetNoAdj(net)
  } else {
    neigh=net$neigh
    diag(neigh)=0
    neigh=graph.adjacency(neigh)

    col = rep("skyblue",net$n)
    col[net$start+1] = "darkorange1"

    set.seed(0)
    plot(neigh,...,layout=net$nodes,edge.arrow.size=0,
         vertex.size=10,main=toupper(net$name),vertex.label=NA,
         vertex.color=col)
  }
}



plotCov<-function(net){
  for(i in 1:ncol(net$Xcov)){
    junk=readline(prompt="Press [enter] to continue\n")
    dat = data.frame(x=net$nodes[,1],
                     y=net$nodes[,2],
                     v=net$Xcov[,i])
    p = ggplot(data=dat,aes(x=x,y=y,fill=v))
    p = p + geom_point(color="black",pch=21,size=10)
    p = p + scale_fill_gradient2(paste(i))
    print(p)
  }
}


ringNetArgs <- function(n){
  bunch = ceiling(n*.05)
  return(list(n-bunch,bunch))
}

ringPlainNetArgs <- function(n){
  return(list(n-1,1))
}

gridNetArgs <- function(n){
  factors = 1:floor(sqrt(n))
  factors = factors[(n %% factors == 0)]
  diff = n
  args = NULL
  for(f0 in factors){
    f1 = n / f0
    if(f1 - f0 < diff){
      diff = f1 - f0
      args = list(f0,f1)
    }
  }
  return(args)
}

randNetArgs <- function(n){
  return(list(n))
}

alleyNetArgs <- function(n){
  prongs = 1
  nextProng = 1
  first = 0
  while(length(prongs) + sum(prongs) <= n){
    if(first == 0){
      prongs = c(prongs,nextProng)
      nextProng = nextProng + 1
      first = 1
    }
    else{
      prongs = c(prongs,nextProng)
      first = 0
    }
  }

  prongs = prongs[-length(prongs)]

  left = n - (length(prongs) + sum(prongs))
  if(left > 0){
    for(i in 0:(left-1)){
      prongs[length(prongs) - i] = prongs[length(prongs) - i] + 1
    }
  }

  args = list()
  for(p in prongs){
    args = c(args,list(rep(1,p)))
  }

  return(list(args))
}


bowTieNetArgs <- function(n){
  randN = ceiling(n*0.1)
  if((n - randN) %% 2 != 0)
    randN = randN + 1
  gridN = (n - randN) / 2

  rand = randNetArgs(randN)
  grid = gridNetArgs(gridN)
  return(list(unlist(grid),unlist(grid),unlist(rand)))
}


scaleFreeNetArgs <- function(n){
  return(list(n))
}


generateNets <- function(n,display=TRUE){
  nets = c("alleyNet","bowTieNet","gridNet","randNet",
           "scaleFreeNet","crpNet")
  nets = foreach(net = nets)%do%{
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


crpNetArgs <- function(n) {
  return(list(n))
}


generateAndSaveNets <- function(nVals){
  nets = c("alleyNet","bowTieNet","gridNet","randNet",
           "scaleFreeNet","crpNet")
  foreach(n = nVals)%:%
    foreach(net = nets)%do%{
      argGen = get(paste(net,"Args",sep=""))
      netGen = get(paste("gen",
                         paste(toupper(substring(net,1,1)),
                               substring(net,2),sep=""),sep=""))

      netRes = do.call(netGen,argGen(n))

      saveNet(netRes)

      cat(paste(net," of size ", n, " is done","\n"))
    }
}


if(FALSE){

  nets100 = generateNets(100,FALSE)
  for(i in nets100)
    saveNet(i)

  nets500 = generateNets(500,FALSE)
  for(i in nets500)
    saveNet(i)

  nets1000 = generateNets(1000,FALSE)
  for(i in nets1000)
    saveNet(i)

  nets10000 = generateNets(10000,FALSE)
  for(i in nets10000)
    saveNet(i)



  grid500 = do.call(genGridNet,args=gridNetArgs(500))
  saveNet(grid500)



  ## shrink net
  shrink=genShrinkNet(17,6)
  saveNet(shrink)



  ## for writeup
  net50 = generateNets(50)
  for(i in net50){
    pdf(paste(i$name,"Sample",".pdf",sep=""))
    plotNet(i,main="")
    dev.off()
  }

  ## for testing
  gridNet25 = do.call(genGridNet,args=gridNetArgs(25))
  saveNet(gridNet25)


  plotN = c(10,100,1000)

  gridPlot = list()
  for(n in plotN){
    cat(paste(n,"\n"))
    gridPlot = c(gridPlot,list(do.call(genGridNet,
                                       args=gridNetArgs(n))))
  }
  for(i in 1:length(plotN)){
    pdf(paste("gridPlot",plotN[i],".pdf",sep=""))
    plotNet(gridPlot[[i]],
            main=paste(plotN[i],"Locations"),
            vertex.size=2)
    dev.off()
  }

  scaleFreePlot = list()
  for(n in plotN){
    cat(paste(n,"\n"))
    scaleFreePlot = c(scaleFreePlot,list(do.call(genScaleFreeNet,
                                                 args=scaleFreeNetArgs(n))))
  }
  for(i in 1:length(plotN)){
    pdf(paste("scaleFreePlot",plotN[i],".pdf",sep=""))
    plotNet(scaleFreePlot[[i]],
            main=paste(plotN[i],"Locations"),
            vertex.size=2)
    dev.off()
  }

  randPlot = list()
  for(n in plotN){
    cat(paste(n,"\n"))
    randPlot = c(randPlot,list(do.call(genRandNet,
                                       args=randNetArgs(n))))
  }
  for(i in 1:length(plotN)){
    pdf(paste("randPlot",plotN[i],".pdf",sep=""))
    plotNet(randPlot[[i]],
            main=paste(plotN[i],"Locations"),
            vertex.size=2)
    dev.off()
  }

}


extract<-function(mat,ind){
  if(ind > 0){
    mat = mat[1:(nrow(mat) - ind + 1),ind:ncol(mat)]
    return(diag(as.matrix(mat)))
  }
  else if(ind < 0){
    ind = abs(ind)
    ind = ind + 1
    mat = mat[1:(nrow(mat) - ind + 1),ind:ncol(mat)]
    return(mat[upper.tri(mat,diag=TRUE)])
  }
}




if(FALSE){
  ## testing the merge cluster routines

  sourceCpp("mergeClusters.cpp")
  cols = colorRampPalette(brewer.pal(9,"Spectral"))

  net = genRandNet(100)

  clusters = getClusters(net$neigh,net$n)
  nClusters = length(unique(clusters))

  plotNet(net,vertex.color=cols(nClusters)[clusters+1])

  while(nClusters > 1){

    net$neigh = mergeOneCluster(net$neigh,clusters,net$nodes[,1],
                              net$nodes[,2],net$n)
    net$neigh = matrix(net$neigh,ncol=net$n)

    readline("Merge clusters")
    plotNet(net,vertex.color=cols(nClusters)[clusters+1])

    clusters = getClusters(net$neigh,net$n)
    nClusters = length(unique(clusters))

    readline("Remove cluster")
    plotNet(net,vertex.color=cols(nClusters)[clusters+1])
  }

}
