#include "rankAgentToyOld.hpp"


template class RankToyOldAgent<GravityModel,GravityParam>;


template <class Model, class ModelParam>
RankToyOldAgent<Model,ModelParam>::RankToyOldAgent(){
  tp.weights.ones(4);
  tp.numChunks = 3;
  tp.valReps=10;
  name="rank";
}

  
template <class Model, class ModelParam>
int RankToyOldAgent<Model,ModelParam>::numFeatures = 4;


template <class Model, class ModelParam>
void RankToyOldAgent<Model,ModelParam>::applyTrt(const SimData & sD,
					   TrtData & tD,
					   const FixedData & fD,
					   const DynamicData & dD,
					   const Model & m,
					   ModelParam & mP){
  if(sD.notInfec.empty())
    return;
  
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);
  
  m.load(sD,tD,fD,dD,mP);
  getSubGraph(fD.numNodes,&fD.network,&subgraph,4);
  
  getFeatures(sD,tD,fD,dD,m,mP);
  std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;
  
  int i,j,node0,addPre,addAct;
  int cI=0,cN=0;
  for(i=0; i<tp.numChunks; i++){

    infRanks = infFeat * tp.weights;
    notRanks = notFeat * tp.weights;


    for(j=0; j<sD.numInfected; j++){
      if(tD.a.at(sD.infected.at(j)))
	sortInfected.push(std::pair<double,int>(std::numeric_limits<double>
						::min(),j));
      else
	sortInfected.push(std::pair<double,int>(infRanks(j),j));
    }
    for(j=0; j<sD.numNotInfec; j++){
      if(tD.p.at(sD.notInfec.at(j)))
	sortNotInfec.push(std::pair<double,int>(std::numeric_limits<double>
						::min(),j));
      else
	sortNotInfec.push(std::pair<double,int>(notRanks(j),j));
    }


    addAct = (int)((i+1)*numAct/std::min(tp.numChunks,numAct)) -
      (int)(i*numAct/std::min(tp.numChunks,numAct));
    for(; cI<(cI+addAct) && cI<numAct; cI++){
      node0=sortInfected.top().second;
      tD.a.at(sD.infected.at(node0)) = 1;
      mP.infProbsBase.row(node0) -= mP.trtAct;
      mP.infProbsSep.row(node0) = 1/(1+arma::exp(mP.infProbsBase.row(node0)));
      sortInfected.pop();
      
    }
    
    addPre = (int)((i+1)*numPre/std::min(tp.numChunks,numPre)) -
      (int)(i*numPre/std::min(tp.numChunks,numPre)); 
    for(; cN<(cN+addPre) && cN<numPre; cN++){
      node0=sortNotInfec.top().second;
      tD.p.at(sD.notInfec.at(node0)) = 1;
      mP.infProbsBase.col(node0) -= mP.trtPre;
      mP.infProbsSep.col(node0) = 1/(1+arma::exp(mP.infProbsBase.col(node0)));
      sortNotInfec.pop();
    }

    
    if((i+1) < tp.numChunks){
      
      getFeatures(sD,tD,fD,dD,m,mP);
    }
  }
}




template <class Model, class ModelParam>
void RankToyOldAgent<Model,ModelParam>::getFeatures(const SimData & sD,
					      const TrtData & tD,
					      const FixedData & fD,
					      const DynamicData & dD,
					      const Model & m,
					      const ModelParam & mP){
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);

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
    
}


std::vector<double> RankToyOldTuneParam::getPar() const {
  std::vector<double> par;
  par = arma::conv_to< std::vector<double> >::from(weights);
  // par.push_back(sigma);
  return par;
}


void RankToyOldTuneParam::putPar(const std::vector<double> & par){
  // sigma = par.back();
  weights = arma::conv_to<arma::colvec>::from(par);
  // weights.resize(weights.n_elem - 1);
}



