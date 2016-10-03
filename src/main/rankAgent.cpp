#include <glog/logging.h>
#include "rankAgent.hpp"

#include <mutex>

std::mutex global_stream_mutex;


template class RankAgent<ToyFeatures5<ModelGravityGDist>,
												 ModelGravityGDist>;

template class RankAgent<ToyFeatures5<Model2GravityGDist>,
												 Model2GravityGDist>;

template class RankAgent<ToyFeatures5<Model2GravityEDist>,
												 Model2GravityEDist>;

template class RankAgent<ToyFeatures5<Model2GPowGDist>,
												 Model2GPowGDist>;

template class RankAgent<ToyFeatures5<Model2EdgeToEdge>,
												 Model2EdgeToEdge>;

template class RankAgent<ToyFeatures5<ModelGDist>,
												 ModelGDist>;

template class RankAgent<ToyFeatures5<ModelEDist>,
												 ModelEDist>;

template class RankAgent<ToyFeatures5<ModelIntercept>,
												 ModelIntercept>;


template class RankAgent<WnsFeatures3<ModelGDist>,
												 ModelGDist>;

template class RankAgent<WnsFeatures3<ModelEDist>,
												 ModelEDist>;

template class RankAgent<WnsFeatures3<ModelIntercept>,
												 ModelIntercept>;

template class RankAgent<WnsFeatures3<Model2EdgeToEdge>,
												 Model2EdgeToEdge>;

template class RankAgent<WnsFeatures3<ModelGravityGDist>,
												 ModelGravityGDist>;

template class RankAgent<WnsFeatures3<Model2GravityGDist>,
												 Model2GravityGDist>;

template class RankAgent<WnsFeatures3<Model2GravityEDist>,
												 Model2GravityEDist>;




template <class F, class M>
RankAgent<F,M>::RankAgent(){
  tp.weights_r.ones(f.numFeatures);
  tp.weights.ones(f.numFeatures);

  tp.jitterScale = 4.0;

  tp.shuffle = false;

	setEdgeToEdge(false);

  name="rank";
}

template <class F, class M>
void RankAgent<F,M>::setEdgeToEdge(const bool edgeToEdge) {
	this->tp.setEdgeToEdge(edgeToEdge);
	this->f.tp.setEdgeToEdge(edgeToEdge);
}

template <class F, class M>
bool RankAgent<F,M>::getEdgeToEdge() const {
	return this->tp.getEdgeToEdge();
}

template <class F, class M>
void RankAgent<F,M>::reset(){
  tp.weights = tp.weights_r;
}


template <class F, class M>
void RankAgent<F,M>::applyTrt(const SimData & sD,
	TrtData & tD,
	const FixedData & fD,
	const DynamicData & dD,
	M & m){
  if(sD.notInfec.empty())
    return;

  CHECK_EQ(getEdgeToEdge(),f.tp.getEdgeToEdge());

  // number of each type of treatment to give
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);

  // precompute data and get baseline features
  f.preCompData(sD,tD,fD,dD,m);
  f.getFeatures(sD,tD,fD,dD,m);

  // jitter the current weights
  arma::colvec jitter;
  arma::mat featStddev;
  jitter.zeros(f.numFeatures);

  int i=0,j=0,node0=0,addPre=0,addAct=0;
  int cI = 0,cN = 0;

  int numChunks = std::log((double)fD.numNodes) + 1.0;

  numChunks = std::min(std::max(numPre,numAct),numChunks);


  for(i = 0; i < numChunks; i++){

    // get jitter
    featStddev.zeros(0,f.numFeatures);
    for(j = 0; j < sD.numNotInfec; ++j){
      if(tD.p.at(sD.notInfec.at(j)) == 0)
				featStddev.insert_rows(0,f.notFeat.row(j));
    }
    for(j = 0; j < sD.numInfected; ++j){
      if(tD.a.at(sD.infected.at(j)) == 0)
				featStddev.insert_rows(0,f.infFeat.row(j));
    }
    featStddev = arma::cov(featStddev);
    // need to add abs() since values close to zero can sometimes
    // result in small negative numbers due to instability
    jitter = arma::sqrt(arma::abs(featStddev.diag()))*calcJitter();

    {
      bool zeroStddev = false;
      for (int featInd = 0; featInd < f.numFeatures; ++featInd) {
        if(featStddev.diag()(featInd) < 1e-8) {
          zeroStddev = true;
        }
      }

      global_stream_mutex.lock();
      CHECK(!zeroStddev)
        << featStddev.diag() << std::endl
        << "inf feat: " << std::endl
        << f.infFeat << std::endl
        << "not feat: " << std::endl
        << f.notFeat << std::endl;
      global_stream_mutex.unlock();
    }

    for(j = 0; j < f.numFeatures; j++)
      jitter(j) *= njm::rnorm01();

    // calculate ranks
    infRanks = f.infFeat * (tp.weights + jitter);
    notRanks = f.notFeat * (tp.weights + jitter);

    // shuffle the node indices
    std::priority_queue<std::pair<double,int> > shufInfected,shufNotInfec;
    for(j = 0; j < sD.numInfected; ++j)
      shufInfected.push(std::pair<double,int>(njm::runif01(),j));
    for(j = 0; j < sD.numNotInfec; ++j)
      shufNotInfec.push(std::pair<double,int>(njm::runif01(),j));

    // sort the locations by their ranks
    // if treated, get lowest value possible
    std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;

    for(j=0; j<sD.numInfected; j++){
      node0 = shufInfected.top().second;
      shufInfected.pop();
      CHECK(std::isfinite(infRanks(node0)))
        << std::endl
        << "rank: " << infRanks(node0) << std::endl
        << "infFeat: " << f.infFeat.row(j)
        << "weights: " << tp.weights.t()
        << "jitter: " << jitter.t()
        << "featStddev: " << featStddev.diag().t()
        << "calcjitter: " << calcJitter();
      if(tD.a.at(sD.infected.at(node0)))
				sortInfected.push(std::pair<double,int>(std::numeric_limits<double>
						::lowest(),node0));
      else
				sortInfected.push(std::pair<double,int>(infRanks(node0),node0));
    }


    for(j=0; j<sD.numNotInfec; j++){
      node0 = shufNotInfec.top().second;
      shufNotInfec.pop();
      CHECK(std::isfinite(notRanks(node0)))
        << std::endl
        << "rank: " << notRanks(node0) << std::endl
        << "notFeat: " << f.notFeat.row(j)
        << "weights: " << tp.weights.t()
        << "jitter: " << jitter.t()
        << "featStddev: " << featStddev.diag().t()
        << "calcjitter: " << calcJitter();
      if(tD.p.at(sD.notInfec.at(node0)))
				sortNotInfec.push(std::pair<double,int>(std::numeric_limits<double>
						::lowest(),node0));
      else
				sortNotInfec.push(std::pair<double,int>(notRanks(node0),node0));
    }


    // randomly select from the top nodes
    std::priority_queue<std::pair<double,int> > selInfected,selNotInfec;
    for(j = 0; j < (numAct - cI); j++){
      if(tp.shuffle){
				selInfected.push(std::pair<double,int>(njm::runif01(),
						sortInfected.top().second));
      }
      else{
				selInfected.push(std::pair<double,int>(sortInfected.top().first,
						sortInfected.top().second));
      }
      sortInfected.pop();
    }

    for(j = 0; j < (numPre - cN); j++){
      if(tp.shuffle){
				selNotInfec.push(std::pair<double,int>(njm::runif01(),
						sortNotInfec.top().second));
      }
      else{
				selNotInfec.push(std::pair<double,int>(sortNotInfec.top().first,
						sortNotInfec.top().second));
      }
      sortNotInfec.pop();
    }

    // number of locations to add treatment too for this iteration
    addPre = 0;
    if(numPre > 0)
      addPre = (int)((i+1)*numPre/std::min(numChunks,numPre)) -
				(int)(i*numPre/std::min(numChunks,numPre));

    addAct = 0;
    if(numAct > 0)
      addAct = (int)((i+1)*numAct/std::min(numChunks,numAct)) -
				(int)(i*numAct/std::min(numChunks,numAct));

    // add active treatment
    for(j = 0; j < addAct && cI < numAct; cI++,j++){
      node0=selInfected.top().second;
      if(tD.a.at(sD.infected.at(node0)) == 1) {
        std::cout << "inf ranks:" << std::endl
                  << infRanks << std::endl
                  << "***********************" << std::endl;
        CHECK_EQ(tD.a.at(sD.infected.at(node0)),0);
      }
      tD.a.at(sD.infected.at(node0)) = 1;
      selInfected.pop();
    }

    // add preventative treatment
    for(j = 0; j < addPre && cN < numPre; cN++,j++){
      node0=selNotInfec.top().second;
      if(tD.p.at(sD.notInfec.at(node0)) == 1) {
        std::cout << "not ranks:" << std::endl
                  << notRanks << std::endl
                  << "***********************" << std::endl;
        CHECK_EQ(tD.p.at(sD.notInfec.at(node0)),0);
      }
      tD.p.at(sD.notInfec.at(node0)) = 1;
      selNotInfec.pop();
    }

    // if more iterations, update features
    if((i+1) < numChunks){
      f.updateFeatures(sD,tD,fD,dD,m);
    }

  }

#ifndef NJM_NO_DEBUG
  int totPre = 0,totAct = 0;
  // check if valid treatments are given to valid locations
  for(i = 0; i < fD.numNodes; i++){
    if(tD.p.at(i) != 1 && tD.p.at(i) != 0){
      std::cout << "Prevenative treatment not 1 or 0"
                << ": " << tD.p.at(i)
                << std::endl;
      throw(1);
    }
    else if(tD.a.at(i) != 1 && tD.a.at(i) != 0){
      std::cout << "Active treatment not 1 or 0"
                << std::endl;
      throw(1);
    }
    else if(tD.a.at(i) == 1 && sD.status.at(i) < 2){
      std::cout << "Not infected receiving active treatment"
                << std::endl;
      throw(1);
    }
    else if(tD.p.at(i) == 1 && sD.status.at(i) >= 2){
      std::cout << "Infected receiving preventative treament"
                << std::endl;
      throw(1);
    }
    else if(tD.a.at(i) == 1)
      totAct++;
    else if(tD.p.at(i) == 1)
      totPre++;
  }

  // check if total number of treatments are correct
  if(totAct != numAct){
    std::cout << "Not correct amount of active treatments."
              << std::endl
              << "Should be " << numAct << " but is " << totAct << "."
              << std::endl
              << "Number of infected nodes is " << sD.numInfected
              << "(" << sD.infected.size() << ")"
              << std::endl;
    throw(1);
  }
  else if(totPre != numPre){
    std::cout << "Not correct amount of preventative treatments."
              << std::endl
              << "Should be " << numPre << " but is " << totPre << "."
              << std::endl
              << "Number of not infected nodes is " << sD.numNotInfec
              << "(" << sD.notInfec.size() << ")"
              << std::endl;
    throw(1);
  }
#endif
}



std::vector<double> RankTuneParam::getPar() const {
  std::vector<double> par;
  par = arma::conv_to< std::vector<double> >::from(weights);
  // par.push_back(sigma);
return par;
}



void RankTuneParam::putPar(const std::vector<double> & par){
  // sigma = par.back();
  weights = arma::conv_to<arma::colvec>::from(par);
  // weights.resize(weights.n_elem - 1);
}


template <class F, class M>
double RankAgent<F,M>::calcJitter(){
  return tp.jitterScale > 0 ? 1.0/tp.jitterScale : 0.0;
}
