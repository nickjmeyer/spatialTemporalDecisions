#include "paramBeta.hpp"


unsigned int ParamBeta::initParsSize(const FixedData & fD){
  return fD.numCovar; // unsigned int literal
}


void ParamBeta::initInternal(const FixedData & fD){
  parsOld = pars;
  covar = fD.covar;
  covarBeta = std::vector<double>(numNodes,0);
  covarBetaOld = covarBeta;
  numNodes = fD.numNodes;
}


void ParamBeta::updateBefore(){
  parsOld = pars;
  covarBetaOld = covarBeta;
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
			const DynamicData & dD) const {
  int i,j;
  std::vector<double>::iterator it;
  for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
    for(j = 0; j < numNodes; ++j, ++it){ // infected
      *it += covarBeta.at(i);
    }
  }
}


void ParamBeta::updateFill(std::vector<double> & probs,
			   const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD) const {
  int i,j;
  double diff;
  std::vector<double>::iterator it;
  for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
    diff = covarBeta.at(i) - covarBetaOld.at(i);
    for(j = 0; j < numNodes; ++j, ++it){ // infected
      *it += diff;
    }
  }
}

