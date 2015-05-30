#include "toyFeatures6.hpp"

std::vector<double> ToyFeatures6TuneParam::getPar() const{
  return std::vector<double>(0);
}


void ToyFeatures6TuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}





template class ToyFeatures6<ModelGravityGDist>;

template class ToyFeatures6<ModelTimeGDist>;

template class ToyFeatures6<ModelTimeExpCavesGDist>;

template class ToyFeatures6<ModelTimeGDistTrendPow>;

template class ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>;

template class ToyFeatures6<ModelTimeExpCavesEDist>;

template class ToyFeatures6<ModelRadius>;

template class ToyFeatures6<ModelGDist>;

template class ToyFeatures6<ModelGDistKern>;

template class ToyFeatures6<ModelCovar>;



template <class M>
void ToyFeatures6<M>::preCompData(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // pre compute stuff

  // load estimated probabilities of infection
  m.modFill(sD,tD,fD,dD);
  m.setQuick(sD,tD,fD,dD);

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
				 (j,m.oneOnOne(*itD1,*itD0,fD.numNodes)));
	
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
void ToyFeatures6<M>::getFeatures(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);



  arma::mat weightMat(m.getQuick().data(),sD.numInfected,sD.numNotInfec,false);

  
  // start feature construction

  
  int i,featNum=0;
  std::vector<int>::const_iterator itD0,itD1,beg;

  // feature 0
  // probability of infection or infecting
  infFeat.col(featNum) = 1 - arma::prod(weightMat,1);
  notFeat.col(featNum) = 1 - arma::prod(weightMat,0).t();
  

  featNum++;


  tDPre = tD;

  arma::colvec notMx = arma::max(notFeat,0).t();
  arma::colvec notMn = arma::min(notFeat,0).t();
  arma::colvec infMx = arma::max(infFeat,0).t();
  arma::colvec infMn = arma::min(infFeat,0).t();

  for(i = 0; i < numFeatures; ++i){
    if((notMx(i) - notMn(i)) > 1e-15){
      notFeat.col(i) = (notFeat.col(i) - notMn(i))/(notMx(i) - notMn(i));
    }
    if((infMx(i) - infMn(i)) > 1e-15){
      infFeat.col(i) = (infFeat.col(i) - infMn(i))/(infMx(i) - infMn(i));
    }
  }

#ifndef NJM_NO_DEBUG
  if(featNum != numFeatures){
    std::cout << "Error: in getFeatures: featNum != numFeatures"
	      << std::endl;
    throw 1;
  }
#endif
}



template <class M>
void ToyFeatures6<M>::updateFeatures(const SimData & sD,
				     const TrtData & tD,
				     const FixedData & fD,
				     const DynamicData & dD,
				     M & m){
  int i,j,num;
  std::pair<int,int> neighOf;

  m.modFill(sD,tD,fD,dD);
  m.setQuick(sD,tD,fD,dD);

  // update neighbor probs
  for(i = 0; i < sD.numNotInfec; i++){
    num=notNeighOfNum.at(i);
    for(j = 0; j < num; j++){
      neighOf = notNeighOf.at(i).at(j);

      notNeigh.at(neighOf.first).at(neighOf.second).second =
	m.oneOnOne(sD.notInfec.at(i),sD.notInfec.at(neighOf.first),
		   fD.numNodes);
    }
  }

  getFeatures(sD,tD,fD,dD,m);
}


