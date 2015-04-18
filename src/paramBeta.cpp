#include "paramBeta.hpp"


unsigned int ParamBeta::initParsSize(const FixedData & fD){
  return fD.numCovar; // unsigned int literal
}


void ParamBeta::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  covar = fD.covar;
  covarBeta = std::vector<double>(numNodes,0);
}


void ParamBeta::updateBefore(){
}


void ParamBeta::updateAfter(){
  int i,j,k,J = parsSize;
  for(i = 0,k = 0; i < numNodes; ++i){
    covarBeta.at(i) = 0;
    for(j = 0; j < J; ++j, ++k){
      covarBeta.at(i) += covar.at(k) * pars.at(j);
    }
  }
}


void ParamBeta::setFill(std::vector<double> & probs,
			const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  int i,j;
  std::vector<double>::iterator it;
  for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
    for(j = 0; j < numNodes; ++j, ++it){ // infected
      *it += covarBeta.at(i);
    }
  }
}


void ParamBeta::modFill(std::vector<double> & probs,
			const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
}

