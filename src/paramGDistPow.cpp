#include "paramGDistPow.hpp"



unsigned int ParamGDistPow::initParsSize(const FixedData & fD){
  return 2U;
}


std::vector<std::string> ParamGDistPow::initNames(){
  return {"alpha","dPow"};
}

std::vector<bool> ParamGDistPow::initToScale(){
  return {true,false};
}



void ParamGDistPow::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  dist = fD.gDist;
  alphaDist = std::vector<double> (fD.numNodes*fD.numNodes,0.0);
}


void ParamGDistPow::updateBefore(){
}


void ParamGDistPow::updateAfter(){
  double alpha = pars.at(0);
  double dPow = pars.at(1);
  int i,I = numNodes * numNodes;
  for(i = 0; i < I; ++i){
    alphaDist[i] = alpha * std::pow(dist[i],dPow);
  }
}


void ParamGDistPow::setFill(std::vector<double> & probs,
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


void ParamGDistPow::modFill(std::vector<double> & probs,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){
}


std::vector<double> ParamGDistPow::partial(const int notNode,
					   const int infNode,
					   const SimData & sD,
					   const TrtData & tD,
					   const FixedData & fD,
					   const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double d = dist[notNode*fD.numNodes + infNode];
  double dPow = std::pow(d,power);
  double dLog = std::log(d);
  return {d, alpha * dLog * dPow};
}


std::vector<double> ParamGDistPow::partial2(const int notNode,
					    const int infNode,
					    const SimData & sD,
					    const TrtData & tD,
					    const FixedData & fD,
					    const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double d = dist[notNode*fD.numNodes + infNode];
  double dPow = std::pow(d,power);
  double dLog = std::log(d);
  return {0.0 , dLog * dPow, dLog * dPow, alpha * dLog * dLog * dPow};
}
