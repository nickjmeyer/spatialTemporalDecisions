#include "wnsFeatures2.hpp"

std::vector<double> WnsFeatures2TuneParam::getPar() const{
  return std::vector<double>(0);
}


void WnsFeatures2TuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}






template class WnsFeatures2<ModelGravityGDist>;

template class WnsFeatures2<ModelGravityEDist>;

template class WnsFeatures2<ModelTimeGDist>;

template class WnsFeatures2<ModelTimeExpCavesGDist>;

template class WnsFeatures2<ModelTimeExpCavesEDist>;

template class WnsFeatures2<ModelRadius>;

template class WnsFeatures2<ModelGDist>;

template class WnsFeatures2<ModelGDistKern>;

template class WnsFeatures2<ModelTimeExpCavesGDistTrendPowCon>;



template <class M>
void WnsFeatures2<M>::preCompData(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // pre compute stuff

  // load estimated probabilities of infection
  m.modFill(sD,tD,fD,dD);
  m.setQuick(sD,tD,fD,dD);

  // extract  connectiviy for not infected
  int i;
  hpddNotInfec.resize(sD.numNotInfec);
  for(i=0; i<sD.numNotInfec; i++)
    hpddNotInfec(i) = fD.hpdd.at(sD.notInfec.at(i));

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
  int j,node,numNeigh = std::min(int(std::log(fD.numNodes)),
				 sD.numNotInfec);
  beg=sD.notInfec.begin();
  for(i=0,itD0=beg; i<sD.numNotInfec; i++,itD0++){
    // sort by distance
    std::priority_queue<std::pair<double,int> > pqd;
    for(j = 0, itD1 = beg; j < sD.numNotInfec; ++j,++itD1){
      pqd.push(std::pair<double,int>(-fD.gDist.at((*itD0)*fD.numNodes
						  + (*itD1)),
				     j));
    }

    // obtain the closest
    std::vector<int> neigh;
    for(j = 0; j < numNeigh; ++j){
      neigh.push_back(pqd.top().second);
      pqd.pop();
    }
    std::sort(neigh.begin(),neigh.end());

    // enter the fill in the containers
    for(j = 0,itD1=neigh.begin(); j < numNeigh; ++j,++itD1){
      node = sD.notInfec.at(*itD1);
      notNeigh.at(i).push_back(std::pair<int,double>
			       (*itD1,
				m.oneOnOne(node,*itD0,fD.numNodes)));

      notNeighOf.at(*itD1).push_back(std::pair<int,int>(i,notNeighNum[i]));

      ++notNeighNum.at(i);
      ++notNeighOfNum.at(*itD1);
    }
  }
}



template <class M>
void WnsFeatures2<M>::getFeatures(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);

  arma::mat weightMat(m.getQuick().data(),sD.numInfected,sD.numNotInfec,false);

  // start feature construction


  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1,beg;


  // feature 0
  // probability of infection or infecting
  infFeat.col(featNum) = 1 - arma::prod(weightMat,1);
  notFeat.col(featNum) = 1 - arma::prod(weightMat,0).t();


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
      modProb *= 1.0 - 1.0/(1.0 + std::exp(neigh.second));
      modProbTot += modProb;
    }
    notFeat(i,featNum) = modProbTot*notFeat(i,0);
  }

  // arma::mat weightMat(sD.numInfected,sD.numNotInfec);
  // double val;
  // int node0;
  // for(i = 0; i < sD.numInfected; ++i){
  //   node0 = sD.infected.at(i);
  //   for(j = 0; j < sD.numNotInfec; ++j){
  //     val = m.oneOnOne(sD.notInfec.at(j),node0,fD.numNodes);
  //     weightMat(i,j) = 1.0 - 1.0/(1.0 + std::exp(val));
  //   }
  // }

  infFeat.col(featNum) = (1.0 - weightMat) * notFeat.col(featNum);


  featNum++;


  // feature 2
  // weighted half plane data depth
  notFeat.col(featNum) = notFeat.col(0) % hpddNotInfec;

  infFeat.col(featNum) = (1.0 - weightMat) * notFeat.col(0);


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
void WnsFeatures2<M>::updateFeatures(const SimData & sD,
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
