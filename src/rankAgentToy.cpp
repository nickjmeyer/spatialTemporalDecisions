#include "rankAgentToy.hpp"


template class RankToyAgent<ToyFeatures0<GravityModel,GravityParam>,
			    GravityModel,GravityParam>;
template class RankToyAgent<ToyFeatures1<GravityModel,GravityParam>,
			    GravityModel,GravityParam>;
template class RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
			    GravityModel,GravityParam>;

template class RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
			    RangeModel,RangeParam>;

template class RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
			    EbolaModel,EbolaParam>;


template <class F, class M, class MP>
RankToyAgent<F,M,MP>::RankToyAgent(){
  tp.weights.ones(4);
  tp.numChunks = 3;

  tp.jitter = 0.25;
  
  f.tp.valReps = 10;
  
  name="rank";
}

  
template <class F, class M, class MP>
void RankToyAgent<F,M,MP>::applyTrt(const SimData & sD,
				    TrtData & tD,
				    const FixedData & fD,
				    const DynamicData & dD,
				    const M & m,
				    MP & mP){
  if(sD.notInfec.empty())
    return;
  
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);
  
  f.preCompData(sD,tD,fD,dD,m,mP);
  f.getFeatures(sD,tD,fD,dD,m,mP);

  
  std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;

  arma::colvec jitter;
  jitter.zeros(f.numFeatures);
  
  int i,j,node0,addPre,addAct;
  int cI=0,cN=0;
  for(i=0; i<tp.numChunks; i++){

    for(j=0; j<f.numFeatures; j++)
      jitter(j) = tp.jitter*njm::rnorm01();
    
    infRanks = f.infFeat * (tp.weights + jitter);
    notRanks = f.notFeat * (tp.weights + jitter);

    sortInfected = std::priority_queue<std::pair<double,int> > ();
    sortNotInfec = std::priority_queue<std::pair<double,int> > ();

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


    addAct = (int)((i+1)*numAct/std::min(tp.numChunks,numAct)) -
      (int)(i*numAct/std::min(tp.numChunks,numAct));
    for(j = 0; j < addAct && cI < numAct; cI++,j++){
      node0=sortInfected.top().second;
      tD.a.at(sD.infected.at(node0)) = 1;
      sortInfected.pop();
      
    }
    
    addPre = (int)((i+1)*numPre/std::min(tp.numChunks,numPre)) -
	     (int)(i*numPre/std::min(tp.numChunks,numPre));
    for(j = 0; j < addPre && cN<numPre; cN++,j++){
      node0=sortNotInfec.top().second;
      tD.p.at(sD.notInfec.at(node0)) = 1;
      sortNotInfec.pop();
    }

    
    if((i+1) < tp.numChunks){
      
      f.updateFeatures(sD,tD,fD,dD,m,mP);
    }
  }
}



std::vector<double> RankToyTuneParam::getPar() const {
  std::vector<double> par;
  par = arma::conv_to< std::vector<double> >::from(weights);
  // par.push_back(sigma);
  return par;
}



void RankToyTuneParam::putPar(const std::vector<double> & par){
  // sigma = par.back();
  weights = arma::conv_to<arma::colvec>::from(par);
  // weights.resize(weights.n_elem - 1);
}



