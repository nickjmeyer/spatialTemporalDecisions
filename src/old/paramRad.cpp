#include "paramRad.hpp"



unsigned int ParamRad::initParsSize(const FixedData & fD){
  return 1U;
}


std::vector<std::string> ParamRad::initNames(){
  return {"rad"};
}


void ParamRad::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  int i,j;
  radVal.clear();
  radVal.reserve(numNodes*numNodes);
  for(i = 0; i < numNodes; ++i){
    for(j = 0; j < numNodes; ++j){
      double n_i = fD.caves[i];
      double n_j = fD.caves[j];
      double cm_ij = fD.cm[i*numNodes + j];
      radVal.push_back(n_i*n_i*n_j/(n_i*cm_ij*(n_i + n_j + cm_ij)));
    }
  }

  strength = std::vector<double>(numNodes*numNodes,0.0);
}


void ParamRad::updateBefore(){
}


void ParamRad::updateAfter(){
  double rad = *beg;
  int i,I = numNodes * numNodes;
  for(i = 0; i < I; ++i){
    strength[i] = rad * radVal[i];
  }
}


void ParamRad::setFill(std::vector<double> & probs,
			 const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD){
  int i,j;
  std::vector<double>::iterator it0;
  for(i = 0,it0 = probs.begin(); i < numNodes; ++i){
    for(j = 0; j < numNodes; ++j,++it0){
      *it0 += strength[i*numNodes + j];
    }
  }
}


void ParamRad::modFill(std::vector<double> & probs,
		       const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD){
}


std::vector<double> ParamRad::partial(const int notNode,
				      const int infNode,
				      const SimData & sD,
				      const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD){
  return std::vector<double>(1,strength[notNode*numNodes+infNode]);
}
