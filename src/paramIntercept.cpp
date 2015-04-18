#include "paramIntercept.hpp"


unsigned int ParamIntercept::initParsSize(const FixedData & fD){
  return 1U; // unsigned int literal
}


void ParamIntercept::initInternal(const FixedData & fD){
}


void ParamIntercept::updateBefore(){
}


void ParamIntercept::updateAfter(){
}


void ParamIntercept::setFill(std::vector<double> & probs,
			     const SimData & sD,
			     const TrtData & tD,
			     const FixedData & fD,
			     const DynamicData & dD){
  double val = *beg;
  std::for_each(probs.begin(),probs.end(),
		[&val](double & x){
		  x += val;
		});
}


void ParamIntercept::modFill(std::vector<double> & probs,
			     const SimData & sD,
			     const TrtData & tD,
			     const FixedData & fD,
			     const DynamicData & dD){
  // do nothing
}
							   
