#include "paramDistKern.hpp"



unsigned int ParamDistKern::initParsSize(const FixedData & fD){
  return 2U;
}


void ParamDistKern::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  dist = fD.dist;
  distKern = std::vector<double> (fD.numNodes*fD.numNodes,0.0);
}


void ParamDistKern::updateBefore(){
}


void ParamDistKern::updateAfter(){
  double alpha = pars.at(0);
  double sigma = pars.at(1);
  int i,I = numNodes * numNodes;
  for(i = 0; i < I; ++i){
    distKern[i] = alpha * std::exp(-dist[i]*dist[i]/(2*std::exp(sigma)));
  }
}


void ParamDistKern::setFill(std::vector<double> & probs,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){
  int i,I = numNodes*numNodes;
  std::vector<double>::iterator it0;
  std::vector<double>::const_iterator it1;
  for(i = 0,
	it0 = probs.begin(),
	it1 = distKern.begin(); i < I; ++it0,++it1,++i){ // not infected
    *it0 -= *it1;
  }
}


void ParamDistKern::modFill(std::vector<double> & probs,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){
}


std::vector<double> ParamDistKern::partial(const int notNode,
					   const int infNode,
					   const SimData & sD,
					   const TrtData & tD,
					   const FixedData & fD,
					   const DynamicData & dD){
  double alpha = pars.at(0);
  double sigma = pars.at(1);
  
  std::vector<double> p2;
  double d = dist[notNode*fD.numNodes + infNode];
  double v = std::exp(-d*d/(2*std::exp(sigma)));
  p2.push_back(-v);
  p2.push_back(-alpha*(d*d/2)*v/std::exp(sigma));
  return p2;
}




std::vector<double> ParamDistKern::partial2(const int notNode,
					    const int infNode,
					    const SimData & sD,
					    const TrtData & tD,
					    const FixedData & fD,
					    const DynamicData & dD){
  double alpha = pars.at(0);
  double sigma = pars.at(1);
  
  std::vector<double> p2;
  double d = dist[notNode*fD.numNodes + infNode];
  double v = std::exp(-d*d/(2*std::exp(sigma)));
  p2.push_back(0);
  p2.push_back(-(d*d/2)*v/std::exp(sigma));
  p2.push_back(-(d*d/2)*v/std::exp(sigma));
  p2.push_back(alpha*(d*d/2)*v/std::exp(sigma)
	       - alpha*(d*d*d*d/4)*v/std::exp(2*sigma));
  return p2;
}
