#include "paramDist.hpp"



unsigned int ParamDist::initParsSize(const FixedData & fD){
  return 1U;
}


void ParamDist::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  dist = fD.dist;
  alphaDist = std::vector<double> (fD.numNodes*fD.numNodes,0.0);
}


void ParamDist::updateBefore(){
}


void ParamDist::updateAfter(){
  double alpha = *beg;
  int i,I = numNodes * numNodes;
  for(i = 0; i < I; ++i){
    alphaDist[i] = alpha * dist[i];
  }
}


void ParamDist::setFill(std::vector<double> & probs,
			const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  int i,I = numNodes*numNodes;
  std::vector<double>::iterator it0;
  std::vector<double>::const_iterator it1;
  for(i = 0,
	it0 = probs.begin(),
	it1 = alphaDist.begin(); i < I; ++it0,++it1,++i){ // not infected
    *it0 -= *it1;
  }
}


void ParamDist::modFill(std::vector<double> & probs,
			const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
}


std::vector<double> ParamDist::partial(const int notNode,
				       const int infNode,
				       const SimData & sD,
				       const TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD){
  return std::vector<double>(1,-dist[notNode*fD.numNodes + infNode]);
}
