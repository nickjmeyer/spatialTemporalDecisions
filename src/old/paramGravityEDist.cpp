#include "paramGravityEDist.hpp"



void ParamGravityEDist::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  grav = std::vector<double>(numNodes*numNodes,0.0);
  
  dist = fD.eDist;
  
  cc.clear();
  cc.reserve(numNodes*numNodes);
  int i,j;
  for(i = 0; i < numNodes; ++i){
    for(j = 0; j < numNodes; ++j){
      cc.push_back(fD.caves.at(i)*fD.caves.at(j));
    }
  }
}
