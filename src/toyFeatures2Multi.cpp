#include "toyFeatures2Multi.hpp"

std::vector<double> ToyFeatures2MultiTuneParam::getPar() const{
  return std::vector<double>(0);
}


void ToyFeatures2MultiTuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}



template class ToyFeatures2Multi<MultiModel>;



template <class M>
int ToyFeatures2Multi<M>::numFeatures = M::numModels*3 + 1;


template <class M>
void ToyFeatures2Multi<M>::preCompData(const SimData & sD,
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
  int S = m.size();
  notNeigh.resize(S);
  std::fill(notNeigh.begin(),notNeigh.end(),
	    std::vector<std::vector<std::pair<int,double> > >(0));
  notNeighOf.resize(sD.numNotInfec);
  std::fill(notNeighOf.begin(),notNeighOf.end(),
	    std::vector<std::pair<int,int> >(0));

  notNeighNum.resize(sD.numNotInfec);
  notNeighOfNum.resize(sD.numNotInfec);
  std::fill(notNeighNum.begin(),notNeighNum.end(),0);
  std::fill(notNeighOfNum.begin(),notNeighOfNum.end(),0);
  
  
  std::vector<int>::const_iterator itD0,itD1,beg;
  int j,s;
  beg=sD.notInfec.begin();
  for(s = 0; s < S; ++s){
    m.modSel(s);

    notNeigh.at(s).resize(sD.numNotInfec);
    std::fill(notNeigh.at(s).begin(),notNeigh.at(s).end(),
	      std::vector<std::pair<int,double> >(0));
    
    for(i=0,itD0=beg; i<sD.numNotInfec; i++,itD0++){
      for(j=0,itD1=beg; j<sD.numNotInfec; j++,itD1++){
	
	if(i!=j && fD.network.at((*itD0)*fD.numNodes + (*itD1))){
	  // neighbors of i
	  notNeigh.at(s).at(i).push_back(std::pair<int,double>
					 (j,m.oneOnOne(*itD1,*itD0,
						       sD,tD,fD,dD)));
	
	  if(s == 0){
	    // i is a neighbor of j
	    notNeighOf.at(j).push_back(std::pair<int,int>(i,notNeighNum.at(i)));

	    // increment totals
	    notNeighNum.at(i)++;
	    notNeighOfNum.at(j)++;
	  }
	  
	}
	
      }
    }
  }
}



template <class M>
void ToyFeatures2Multi<M>::getFeatures(const SimData & sD,
				       const TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD,
				       M & m){
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);


  
  // start feature construction

  int s,S = m.size();
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1,beg;

  arma::colvec infProbs;
  
  for(s = 0; s < S; ++s){ // begin model dependent features
    m.modSel(s);
    
    // feature 0
    // probability of infection or infecting
    infFeat.col(featNum) = 1 - arma::prod(m.getPar()->getSep(),1);
    notFeat.col(featNum) = 1 - arma::prod(m.getPar()->getSep(),0).t();

    infProbs = notFeat.col(featNum);
  
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
	neigh=notNeigh.at(s).at(i).at(j);
      
	modProb = 1.0 - infProbs(neigh.first);
	modProb *= 1.0/(1.0 + std::exp(neigh.second));
	modProbTot += 1.0 - modProb;
      }
      notFeat(i,featNum) = modProbTot*infProbs(i);
    }
  
    infFeat.col(featNum) = (1.0-m.getPar()->getSep()) * notFeat.col(featNum);

    featNum++;

  
    // feature 2
    // weighted subgraph connectivity measures
    notFeat.col(featNum) = infProbs % subGraphNotInfec;

    infFeat.col(featNum) = (1.0 - m.getPar()->getSep()) * infProbs;

    featNum++;

  } // end model dependent features
  

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

#ifndef NJM_NO_DEBUG
  if(featNum != numFeatures){
    std::cout << "Error: in getFeatures: featNum != numFeatures"
	      << std::endl;
    throw 1;
  }
#endif
    
}



template <class M>
void ToyFeatures2Multi<M>::updateFeatures(const SimData & sD,
					  const TrtData & tD,
					  const FixedData & fD,
					  const DynamicData & dD,
					  M & m){
  int i,j,node0,now,pre,num;
  std::pair<int,int> neighOf;

  int s,S = m.size();

  for(s = 0; s < S; ++s){
    m.modSel(s);
  
    // update not infected probabilities
    for(i = 0; i < sD.numNotInfec; i++){
      node0 = sD.notInfec.at(i);
      now = tD.p.at(node0);
      pre = tDPre.p.at(node0);
    
      if(now != pre && now == 1){ // adding trt
	m.getPar()->getBase().col(i) -= m.getPar()->getTrtPre();
	m.getPar()->setCol(i);

	num=notNeighOfNum.at(i);
	for(j = 0; j < num; j++){
	  neighOf = notNeighOf.at(i).at(j);

	  notNeigh.at(s).at(neighOf.first).at(neighOf.second).second
	    -= m.getPar()->getTrtPre();
	}
      }
    
      else if(now != pre && now == 0){ // removing trt
	m.getPar()->getBase().col(i) += m.getPar()->getTrtPre();
	m.getPar()->setCol(i);

	num=notNeighOfNum.at(i);
	for(j = 0; j < num; j++){
	  neighOf = notNeighOf.at(i).at(j);

	  notNeigh.at(s).at(neighOf.first).at(neighOf.second).second
	    += m.getPar()->getTrtPre();
	}
      }
    }

    
    // update infected probabilities
    for(i = 0; i < sD.numInfected; i++){
      node0 = sD.infected.at(i);
      now = tD.a.at(node0);
      pre = tDPre.a.at(node0);
    
      if(now != pre && now == 1){ // adding trt
	m.getPar()->getBase().row(i) -= m.getPar()->getTrtAct();
	m.getPar()->setRow(i);
      }
      else if(now != pre && now == 0){ // removing trt
	m.getPar()->getBase().row(i) += m.getPar()->getTrtAct();
	m.getPar()->setRow(i);
      }
    }

  }

  getFeatures(sD,tD,fD,dD,m);
}


