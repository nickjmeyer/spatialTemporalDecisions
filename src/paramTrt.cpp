#include "paramTrt.hpp"



unsigned int ParamTrt::initParsSize(const FixedData & fD){
  return 2U;
}


void ParamTrt::initInternal(const FixedData & fD){
  parsOld = pars;
  numNodes = fD.numNodes;

  a = std::vector<int>(numNodes,0);
  p = std::vector<int>(numNodes,0);
  trt.clear();
}


void ParamTrt::updateBefore(){
  parsOld = pars;
  trtOld = trt;
}


void ParamTrt::updateAfter(){
  double trtAct = pars.at(0);
  double trtPre = pars.at(1);
  int i,j,k,I = trt.size();
  for(i = 0; i < I; ++i){
    j = trt.at(i).first / numNodes; // not infected
    k = trt.at(i).first % numNodes; // infected
    trt.at(i).second = trtAct * double(a.at(k)) + trtPre * double(p.at(k));
  }
}


void ParamTrt::setFill(std::vector<double> & probs,
		       const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD) const {
  int i,I = trt.size();
  std::vector<double>::iterator it0;
  std::vector<std::pair<int,double> >::const_iterator it1;
  for(i = 0,
	it0 = probs.begin(),
	it1 = trt.begin(); i < I; ++it0,++it1,++i){ // not infected
    *it0 += (*it1).second;
  }
}


void ParamTrt::updateFill(std::vector<double> & probs,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const {
  int i,I = trt.size();
  std::vector<double>::iterator it0;
  std::vector<std::pair<int,double> >::const_iterator it1,it2;
  for(i = 0,
	it0 = probs.begin(),
	it1 = trt.begin(); i < I; ++it0,++it1,++i){ // not infected
    *it0 += (*it1).second;
  }
  for(i = 0, it0 = probs.begin(),
	it1 = trt.begin(),
	it2 = trtOld.begin(); i < I; ++it0,++it1,++it2,++i){ // not infected
    *it0 += (*it1).second - (*it2).second;
  }
}

