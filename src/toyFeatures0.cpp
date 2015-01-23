#include "toyFeatures0.hpp"

std::vector<double> ToyFeatures0TuneParam::getPar() const{
  return std::vector<double>(0);
}


void ToyFeatures0TuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}





template class ToyFeatures0<GravityModel,GravityParam>;




template <class Model, class ModelParam>
int ToyFeatures0<Model,ModelParam>::numFeatures = 4;


template <class Model, class ModelParam>
void ToyFeatures0<Model,ModelParam>::preCompData(const SimData & sD,
						 const TrtData & tD,
						 const FixedData & fD,
						 const DynamicData & dD,
						 const Model & m,
						 ModelParam & mP){
  // pre compute stuff
  m.load(sD,tD,fD,dD,mP);
  getSubGraph(fD.numNodes,&fD.network,&subgraph,4);
}



template <class Model, class ModelParam>
void ToyFeatures0<Model,ModelParam>::getFeatures(const SimData & sD,
						 const TrtData & tD,
						 const FixedData & fD,
						 const DynamicData & dD,
						 const Model & m,
						 ModelParam & mP){
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);


  
  // start feature construction


  
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;


  // feature 0
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,0).t();
  
  featNum++;
  
  // feature 1
  SystemLight<Model,ModelParam> s(sD,tD,fD,dD,m,mP);
  std::vector<int> newInfec;
  int k,numNewInfec;
  for(i=0; i<tp.valReps; i++){
    s.reset();
    s.nextPoint();
    
    newInfec=s.sD.newInfec;
    s.nextPoint(1);
    
    newInfec.insert(newInfec.end(),s.sD.newInfec.begin(),s.sD.newInfec.end());
    
    numNewInfec = newInfec.size();
    for(j=0,itD0=newInfec.begin(); j<numNewInfec; j++,itD0++)
      for(k=0,itD1=sD.notInfec.begin(); k<sD.numNotInfec; k++,itD1++)
	if(*itD0 == *itD1)
	  notFeat(k,featNum)+=1.0/(double)tp.valReps;
  }

  infFeat.col(featNum) = (1.0-mP.infProbsSep) * notFeat.col(featNum);
  
  featNum++;

  
  // feature 2
  std::vector<double> subgraphNoTrt;

  std::vector<int> networkNoTrt;
  int noneTrt=1;
  i=0;
  while(noneTrt && i < fD.numNodes){
    if(tD.a.at(i) || tD.p.at(i))
      noneTrt=0;
    i++;
  }
    
  if(!noneTrt){
    networkNoTrt=fD.network;
    for(i=0; i<fD.numNodes; i++){
      if(tD.a.at(i) || tD.p.at(i))
	for(j=0; j<fD.numNodes; j++){
	  networkNoTrt.at(i*fD.numNodes + j) = 0;
	  networkNoTrt.at(j*fD.numNodes + i) = 0;
	}
    }
    getSubGraph(fD.numNodes,&networkNoTrt,&subgraphNoTrt,4);
  }

  int node0;
  double subgraphNode,minDist;
  if(!noneTrt){
    for(i=0; i<sD.numNotInfec; i++){
      node0=sD.notInfec.at(i);
      minDist=0;
      for(j=0; j<sD.numInfected; j++)
	if(1.0/fD.logDist.at(node0*fD.numNodes + sD.infected.at(j)) > minDist)
	  minDist = 1.0/fD.logDist.at(node0*fD.numNodes + sD.infected.at(j));
      subgraphNode=subgraph.at(node0);
      if(subgraphNode>0)
	notFeat(i,featNum) = minDist*(subgraphNode - subgraphNoTrt.at(node0))/
	  (subgraphNode);
    }
    for(i=0; i<sD.numInfected; i++){
      node0=sD.infected.at(i);
      minDist=0;
      for(j=0; j<sD.numNotInfec; j++)
	if(1.0/fD.logDist.at(node0*fD.numNodes + sD.notInfec.at(j)) > minDist)
	  minDist = 1.0/fD.logDist.at(node0*fD.numNodes + sD.notInfec.at(j));
      subgraphNode=subgraph.at(node0);
      if(subgraphNode>0)
	infFeat(i,featNum) = minDist*(subgraphNode - subgraphNoTrt.at(node0))/
	  (subgraphNode);
    }

  }
  else{
    double maxSubgraph = 0;
    for(i=0; i<fD.numNodes; i++)
      if(subgraph.at(i) > maxSubgraph)
	maxSubgraph = subgraph.at(i);
    for(i=0; i<sD.numNotInfec; i++){
      node0=sD.notInfec.at(i);
      minDist=0;
      for(j=0; j<sD.numInfected; j++)
	if(1.0/fD.logDist.at(node0*fD.numNodes + sD.infected.at(j)) > minDist)
	  minDist = 1.0/fD.logDist.at(node0*fD.numNodes + sD.infected.at(j));
      notFeat(i,featNum) = minDist*subgraph.at(node0)/maxSubgraph;
    }
    for(i=0; i<sD.numInfected; i++){
      node0=sD.infected.at(i);
      minDist=0;
      for(j=0; j<sD.numNotInfec; j++)
	if(1.0/fD.logDist.at(node0*fD.numNodes + sD.notInfec.at(j)) > minDist)
	  minDist = 1.0/fD.logDist.at(node0*fD.numNodes + sD.notInfec.at(j));
      infFeat(i,featNum) = minDist*subgraph.at(node0)/maxSubgraph;
    }
  }
    
  featNum++;
  

  // feature 3
  std::vector<int>::const_iterator itD2,itD3;
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    notFeat(i,featNum) = std::log(1.0+totalDist);
    // p=std::priority_queue<double>();
    // for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
    //   p.push(1.0/(1.0+fD.dist.at((*itD0)*fD.numNodes + *itD1)));
    // notFeat(i,featNum) = p.top()/fD.invDistSD;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    infFeat(i,featNum) = std::log(1.0+totalDist);
    // p=std::priority_queue<double>();
    // for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
    //   p.push(1.0/(1.0+fD.dist.at((*itD1)*fD.numNodes + *itD0)));
    // infFeat(i,featNum) = p.top()/fD.invDistSD;
  }

  
  featNum++;

  tDPre = tD;

  if(featNum != numFeatures){
    std::cout << "Error: in getFeatures: featNum != numFeatures"
	      << std::endl;
    throw 1;
  }
    
}



template <class Model, class ModelParam>
void ToyFeatures0<Model,ModelParam>::updateFeatures(const SimData & sD,
						    const TrtData & tD,
						    const FixedData & fD,
						    const DynamicData & dD,
						    const Model & m,
						    ModelParam & mP){
  int i,node0,now,pre;

  // update not infected probabilities
  for(i=0; i<sD.numNotInfec; i++){
    node0 = sD.notInfec.at(i);
    now = tD.p.at(node0);
    pre = tDPre.p.at(node0);
    
    if(now != pre && pre == 1){ // adding trt
      mP.infProbsBase.col(i) -= mP.trtPre;
      mP.infProbsSep.col(i) = 1/(1+arma::exp(mP.infProbsBase.col(i)));
    } 
    else if(now != pre && pre == 0){ // removing trt
      mP.infProbsBase.col(i) += mP.trtPre;
      mP.infProbsSep.col(i) = 1/(1+arma::exp(mP.infProbsBase.col(i)));
    }
  }

  // update infected probabilities
  for(i=0; i<sD.numInfected; i++){
    node0 = sD.infected.at(i);
    now = tD.a.at(node0);
    pre = tDPre.a.at(node0);
    
    if(now != pre && pre == 1){ // adding trt
      mP.infProbsBase.row(i) -= mP.trtAct;
      mP.infProbsSep.row(i) = 1/(1+arma::exp(mP.infProbsBase.row(i)));
    }
    else if(now != pre && pre == 0){ // removing trt
      mP.infProbsBase.row(i) += mP.trtPre;
      mP.infProbsSep.row(i) = 1/(1+arma::exp(mP.infProbsBase.row(i)));
    }
  }


  getFeatures(sD,tD,fD,dD,m,mP);
}


  
