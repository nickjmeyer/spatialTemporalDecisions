#include "paramGravity.hpp"



unsigned int ParamGravity::initParsSize(const FixedData & fD){
  return 2U;
}


void ParamGravity::initInternal(const FixedData & fD){
  parsOld = pars;
  numNodes = fD.numNodes;
  grav = std::vector<double>(numNodes*numNodes,0.0);
  gravOld = grav;
  
  dist = fD.dist;
  
  cc.clear();
  cc.reserve(numNodes*numNodes);
  int i,j;
  for(i = 0; i < numNodes; ++i){
    for(j = 0; j < numNodes; ++j){
      cc.push_back(fD.caves.at(i)*fD.caves.at(j));
    }
  }
}


void ParamGravity::updateBefore(){
  parsOld = pars;
  gravOld = grav;
}


void ParamGravity::updateAfter(){
  double alpha = pars.at(0);
  double power = pars.at(1);
  int i,I = numNodes * numNodes;
  for(i = 0; i < I; ++i){
    grav.at(i) = alpha * dist.at(i) / std::pow(cc.at(i),power);
  }
}


void ParamGravity::setFill(std::vector<double> & probs,
			   const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD) const {
  int i,I = numNodes*numNodes;
  std::vector<double>::iterator it0;
  std::vector<double>::const_iterator it1;
  for(i = 0,
	it0 = probs.begin(),
	it1 = grav.begin(); i < I; ++it0,++it1,++i){ // not infected
    *it0 += *it1;
  }
}


void ParamGravity::updateFill(std::vector<double> & probs,
			      const SimData & sD,
			      const TrtData & tD,
			      const FixedData & fD,
			      const DynamicData & dD) const {
  int i,I = numNodes*numNodes;
  std::vector<double>::iterator it0;
  std::vector<double>::const_iterator it1,it2;
  for(i = 0, it0 = probs.begin(),
	it1 = grav.begin(),
	it2 = gravOld.begin(); i < I; ++it0,++it1,++it2,++i){ // not infected
    *it0 += *it1 - *it2;
  }
}

