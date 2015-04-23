#include "toyFeatures2.hpp"

std::vector<double> ToyFeatures2TuneParam::getPar() const{
  return std::vector<double>(0);
}


void ToyFeatures2TuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}





template class ToyFeatures2<GravityTimeInfExpCavesModel>;

template class ToyFeatures2<GravityTimeInfExpModel>;

template class ToyFeatures2<GravityTimeInfModel>;

template class ToyFeatures2<GravityModel>;

template class ToyFeatures2<RangeModel>;

template class ToyFeatures2<RadiusModel>;

template class ToyFeatures2<CaveModel>;




template <class M>
void ToyFeatures2<M>::preCompData(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // pre compute stuff

  // load estimated probabilities of infection
  m.load(sD,tD,fD,dD);

  // extract subgraph connectiviy for not infected
  int i;
  subGraphNotInfec.resize(sD.numNotInfec);
  for(i=0; i<sD.numNotInfec; i++)
    subGraphNotInfec(i) = fD.subGraph.at(sD.notInfec.at(i));

  // obtain neighbors and probabilities not infected infects other not infected
  // initialize containers
  notNeigh.resize(sD.numNotInfec);
  notNeighOf.resize(sD.numNotInfec);
  std::fill(notNeigh.begin(),notNeigh.end(),
	    std::vector<std::pair<int,double> >(0));
  std::fill(notNeighOf.begin(),notNeighOf.end(),
	    std::vector<std::pair<int,int> >(0));

  notNeighNum.resize(sD.numNotInfec);
  notNeighOfNum.resize(sD.numNotInfec);
  std::fill(notNeighNum.begin(),notNeighNum.end(),0);
  std::fill(notNeighOfNum.begin(),notNeighOfNum.end(),0);
  
  
  std::vector<int>::const_iterator itD0,itD1,beg;
  int j;
  beg=sD.notInfec.begin();
  for(i=0,itD0=beg; i<sD.numNotInfec; i++,itD0++){
    for(j=0,itD1=beg; j<sD.numNotInfec; j++,itD1++){
      if(i!=j && fD.network.at((*itD0)*fD.numNodes + (*itD1))){
	// neighbors of i
	notNeigh.at(i).push_back(std::pair<int,double>
				 (j,m.oneOnOne(*itD1,*itD0,sD,tD,fD,dD)));
	
	// i is a neighbor of j
	notNeighOf.at(j).push_back(std::pair<int,int>(i,notNeighNum.at(i)));

	// increment totals
	notNeighNum.at(i)++;
	notNeighOfNum.at(j)++;
      }
    }
  }
}



template <class M>
void ToyFeatures2<M>::getFeatures(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);


  
  // start feature construction

  
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1,beg;

  
  // feature 0
  // probability of infection or infecting
  infFeat.col(featNum) = 1 - arma::prod(m.mP.infProbsSep,1);
  notFeat.col(featNum) = 1 - arma::prod(m.mP.infProbsSep,0).t();
  

  featNum++;

  
  
  // feature 1
  // joint probability of infection between not infected and not infected
  // weighted average of joint probabilities
  int num;
  double modProbTot,modProb;
  std::pair<int,double> neigh;
  for(i = 0; i < sD.numNotInfec; i++){
    modProbTot = 0;
    num = notNeighNum.at(i);
    for(j = 0; j < num; j++){
      neigh=notNeigh.at(i).at(j);
      
      modProb = 1.0 - notFeat(neigh.first,0);
      modProb *= 1.0/(1.0 + std::exp(neigh.second));
      modProbTot += 1.0 - modProb;
    }
    notFeat(i,featNum) = modProbTot*notFeat(i,0);
  }
  
  infFeat.col(featNum) = (1.0-m.mP.infProbsSep) * notFeat.col(featNum);

  
  featNum++;

  
  // feature 2
  // weighted subgraph connectivity measures
  notFeat.col(featNum) = notFeat.col(0) % subGraphNotInfec;

  infFeat.col(featNum) = (1.0 - m.mP.infProbsSep) * notFeat.col(0);

    
  featNum++;


  // feature 3
  // density estimate
  std::vector<int>::const_iterator itD2,itD3;
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  double totalDist;
  for(i = 0,itD0 = itD2; i < sD.numNotInfec; i++, itD0++){
    // density estimates for not infected
    totalDist=0;
    for(j = 0,itD1 = itD3; j < sD.numInfected; j++, itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    notFeat(i,featNum) = std::log(1.0+totalDist);
  }

  for(i = 0,itD0 = itD3; i < sD.numInfected; i++, itD0++){
    // density estimates for infected
    totalDist=0;
    for(j = 0,itD1 = itD2; j < sD.numNotInfec; j++, itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    infFeat(i,featNum) = std::log(1.0+totalDist);
  }

  featNum++;

  tDPre = tD;

  // std::cout << "infProbs: " << arma::sum(infFeat,0)
  // 	    << "notProbs: " << arma::sum(notFeat,0);


#ifndef NJM_NO_DEBUG
  if(featNum != numFeatures){
    std::cout << "Error: in getFeatures: featNum != numFeatures"
	      << std::endl;
    throw 1;
  }
#endif
}



template <class M>
void ToyFeatures2<M>::updateFeatures(const SimData & sD,
				     const TrtData & tD,
				     const FixedData & fD,
				     const DynamicData & dD,
				     M & m){
  int i,j,node0,now,pre,num;
  std::pair<int,int> neighOf;

  // update not infected probabilities
  for(i = 0; i < sD.numNotInfec; i++){
    node0 = sD.notInfec.at(i);
    now = tD.p.at(node0);
    pre = tDPre.p.at(node0);
    
    if(now != pre && now == 1){ // adding trt
      m.mP.infProbsBase.col(i) -= m.mP.trtPre;
      m.mP.setCol(i);

      num=notNeighOfNum.at(i);
      for(j = 0; j < num; j++){
	neighOf = notNeighOf.at(i).at(j);

	notNeigh.at(neighOf.first).at(neighOf.second).second -= m.mP.trtPre;
      }
    }
    
    else if(now != pre && now == 0){ // removing trt
      m.mP.infProbsBase.col(i) += m.mP.trtPre;
      m.mP.setCol(i);

      num=notNeighOfNum.at(i);
      for(j = 0; j < num; j++){
	neighOf = notNeighOf.at(i).at(j);

	notNeigh.at(neighOf.first).at(neighOf.second).second += m.mP.trtPre;
      }
    }
  }

  // update infected probabilities
  for(i = 0; i < sD.numInfected; i++){
    node0 = sD.infected.at(i);
    now = tD.a.at(node0);
    pre = tDPre.a.at(node0);
    
    if(now != pre && now == 1){ // adding trt
      m.mP.infProbsBase.row(i) -= m.mP.trtAct;
      m.mP.setRow(i);
    }
    else if(now != pre && now == 0){ // removing trt
      m.mP.infProbsBase.row(i) += m.mP.trtAct;
      m.mP.setRow(i);
    }
  }
  getFeatures(sD,tD,fD,dD,m);
}


