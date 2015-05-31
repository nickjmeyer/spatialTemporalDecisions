#include "toyFeatures7.hpp"

std::vector<double> ToyFeatures7TuneParam::getPar() const{
  return std::vector<double>(0);
}


void ToyFeatures7TuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}





template class ToyFeatures7<ModelGravityGDist>;

template class ToyFeatures7<ModelTimeGDist>;

template class ToyFeatures7<ModelTimeExpCavesGDist>;

template class ToyFeatures7<ModelTimeGDistTrendPow>;

template class ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>;

template class ToyFeatures7<ModelTimeExpCavesEDist>;

template class ToyFeatures7<ModelRadius>;

template class ToyFeatures7<ModelGDist>;

template class ToyFeatures7<ModelGDistKern>;

template class ToyFeatures7<ModelCovar>;



template <class M>
void ToyFeatures7<M>::preCompData(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
}



template <class M>
void ToyFeatures7<M>::getFeatures(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);



  // start feature construction

  
  int i,j,featNum=0;

  // feature 0
  for(i = 0; i < sD.numInfected; ++i){
    std::priority_queue<double> pq;
    for(j = 0; j < sD.numNotInfec; ++j){
      pq.push(-fD.gDist.at(sD.infected.at(i)*fD.numNodes + sD.notInfec.at(j)));
    }
    infFeat(i,featNum) = -pq.top();
  }

  for(i = 0; i < sD.numNotInfec; ++i){
    std::priority_queue<double> pq;
    for(j = 0; j < sD.numInfected; ++j){
      pq.push(-fD.gDist.at(sD.notInfec.at(i)*fD.numNodes + sD.infected.at(j)));
    }
    notFeat(i,featNum) = -pq.top();
  }
  

  featNum++;


  // feature 1
  for(i = 0; i < sD.numInfected; ++i){
    std::priority_queue<double> pq;
    for(j = 0; j < sD.numNotInfec; ++j){
      pq.push(-fD.gDist.at(sD.infected.at(i)*fD.numNodes + sD.notInfec.at(j)));
    }
    infFeat(i,featNum) = (-pq.top())*(-pq.top());
  }

  for(i = 0; i < sD.numNotInfec; ++i){
    std::priority_queue<double> pq;
    for(j = 0; j < sD.numInfected; ++j){
      pq.push(-fD.gDist.at(sD.notInfec.at(i)*fD.numNodes + sD.infected.at(j)));
    }
    notFeat(i,featNum) = (-pq.top())*(-pq.top());
  }
  

  featNum++;
  

  tDPre = tD;

  // arma::colvec notMx = arma::max(notFeat,0).t();
  // arma::colvec notMn = arma::min(notFeat,0).t();
  // arma::colvec infMx = arma::max(infFeat,0).t();
  // arma::colvec infMn = arma::min(infFeat,0).t();

  // for(i = 0; i < numFeatures; ++i){
  //   if((notMx(i) - notMn(i)) > 1e-15){
  //     notFeat.col(i) = (notFeat.col(i) - notMn(i))/(notMx(i) - notMn(i));
  //   }
  //   if((infMx(i) - infMn(i)) > 1e-15){
  //     infFeat.col(i) = (infFeat.col(i) - infMn(i))/(infMx(i) - infMn(i));
  //   }
  // }

#ifndef NJM_NO_DEBUG
  if(featNum != numFeatures){
    std::cout << "Error: in getFeatures: featNum != numFeatures"
	      << std::endl;
    throw 1;
  }
#endif
}



template <class M>
void ToyFeatures7<M>::updateFeatures(const SimData & sD,
				     const TrtData & tD,
				     const FixedData & fD,
				     const DynamicData & dD,
				     M & m){
}


