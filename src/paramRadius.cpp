#include "paramRadius.hpp"



unsigned int ParamRadius::initParsSize(const FixedData & fD){
  return 1U;
}


void ParamRadius::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  logDist = fD.logDist;
  beyond = std::vector<int> (fD.numNodes*fD.numNodes,0);
}


void ParamRadius::updateBefore(){
}


void ParamRadius::updateAfter(){
  double radius = *beg;
  int i,I = numNodes * numNodes;
  for(i = 0; i < I; ++i){ 
    if(logDist.at(i) > radius)
      beyond.at(i) = 1;
    else
      beyond.at(i) = 0;
  }
}


void ParamRadius::setFill(std::vector<double> & probs,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
  int i,I = numNodes*numNodes;
  std::vector<double>::iterator it0;
  std::vector<int>::const_iterator it1;
  for(i = 0,
	it0 = probs.begin(),
	it1 = beyond.begin(); i < I; ++it0,++it1,++i){ // not infected
    if(*it1 == 1)
      *it0 -= 500.0;
  }
}


void ParamRadius::modFill(std::vector<double> & probs,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
}

