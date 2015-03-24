#include "rankAgent.hpp"


template class RankAgent<ToyFeatures2<GravityTimeInfExpCavesModel,
				      GravityTimeInfExpCavesParam>,
			 GravityTimeInfExpCavesModel,
			 GravityTimeInfExpCavesParam>;

template class RankAgent<ToyFeatures2<GravityTimeInfExpModel,
				      GravityTimeInfExpParam>,
			 GravityTimeInfExpModel,
			 GravityTimeInfExpParam>;

template class RankAgent<ToyFeatures2<GravityTimeInfModel,GravityTimeInfParam>,
			 GravityTimeInfModel,GravityTimeInfParam>;

template class RankAgent<ToyFeatures2<GravityModel,GravityParam>,
			 GravityModel,GravityParam>;

template class RankAgent<ToyFeatures2<RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class RankAgent<ToyFeatures2<RadiusModel,RadiusParam>,
			 RadiusModel,RadiusParam>;

template class RankAgent<ToyFeatures2<CaveModel,CaveParam>,
			 CaveModel,CaveParam>;



template <class F, class M, class MP>
RankAgent<F,M,MP>::RankAgent(){
  tp.weights_r.ones(f.numFeatures);
  tp.weights.ones(f.numFeatures);

  tp.jitterScale = 4.0;
  
  name="rank";
}


template <class F, class M, class MP>
void RankAgent<F,M,MP>::reset(){
  tp.weights = tp.weights_r;
}

  
template <class F, class M, class MP>
void RankAgent<F,M,MP>::applyTrt(const SimData & sD,
				 TrtData & tD,
				 const FixedData & fD,
				 const DynamicData & dD,
				 const M & m,
				 MP & mP){
  if(sD.notInfec.empty())
    return;

  // number of each type of treatment to give
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);

  // precompute data and get baseline features
  f.preCompData(sD,tD,fD,dD,m,mP);
  f.getFeatures(sD,tD,fD,dD,m,mP);

  // jitter the current weights
  arma::colvec jitter;
  arma::mat featStddev;
  jitter.zeros(f.numFeatures);
  
  int i,j,node0,addPre,addAct;
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
    jitter = arma::sqrt(featStddev.diag())/tp.jitterScale;
    
    for(j = 0; j < f.numFeatures; j++)
      jitter(j) *= njm::rnorm01();

    // calculate ranks
    infRanks = f.infFeat * (tp.weights + jitter);
    notRanks = f.notFeat * (tp.weights + jitter);


    // sort the locations by their ranks
    // if treated, get lowest value possible
    std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;
    
    for(j=0; j<sD.numInfected; j++){
      if(tD.a.at(sD.infected.at(j)))
	sortInfected.push(std::pair<double,int>(std::numeric_limits<double>
						::lowest(),j));
      else
	sortInfected.push(std::pair<double,int>(infRanks(j),j));
    }
    
    for(j=0; j<sD.numNotInfec; j++){
      if(tD.p.at(sD.notInfec.at(j)))
	sortNotInfec.push(std::pair<double,int>(std::numeric_limits<double>
						::lowest(),j));
      else
	sortNotInfec.push(std::pair<double,int>(notRanks(j),j));
    }


    std::priority_queue<std::pair<double,int> > selInfected,selNotInfec;
    for(j = 0; j < (numAct - cI); j++){
      selInfected.push(std::pair<double,int>(njm::runif01(),
					     sortInfected.top().second));
      sortInfected.pop();
    }

    for(j = 0; j < (numPre - cN); j++){
      selNotInfec.push(std::pair<double,int>(njm::runif01(),
					     sortNotInfec.top().second));
      sortNotInfec.pop();
    }
    

    // number of locations to add treatment too for this iteration
    addPre = (int)((i+1)*numPre/std::min(numChunks,numPre)) -
      (int)(i*numPre/std::min(numChunks,numPre));
    addAct = (int)((i+1)*numAct/std::min(numChunks,numAct)) -
      (int)(i*numAct/std::min(numChunks,numAct));


    // add active treatment
    for(j = 0; j < addAct && cI < numAct; cI++,j++){
      node0=selInfected.top().second;
      tD.a.at(sD.infected.at(node0)) = 1;
      selInfected.pop();
    }

    // add preventative treatment
    for(j = 0; j < addPre && cN < numPre; cN++,j++){
      node0=selNotInfec.top().second;
      tD.p.at(sD.notInfec.at(node0)) = 1;
      selNotInfec.pop();
    }

    // if more iterations, update features
    if((i+1) < numChunks){
      f.updateFeatures(sD,tD,fD,dD,m,mP);
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



