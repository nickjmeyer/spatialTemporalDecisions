#include "paramTrendPow.hpp"


unsigned int ParamTrendPow::initParsSize(const FixedData & fD){
  return 2U; // unsigned int literal
}


void ParamTrendPow::initInternal(const FixedData & fD){
  prevTime = 0;
}


void ParamTrendPow::updateBefore(){
}


void ParamTrendPow::updateAfter(){
}


void ParamTrendPow::setFill(std::vector<double> & probs,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double val = alpha*std::pow(double(sD.time),-std::exp(power));
  prevTime = (unsigned int)(sD.time);
  std::for_each(probs.begin(),probs.end(),
		[&val](double & x){
		  x += val;
		});
}


void ParamTrendPow::modFill(std::vector<double> & probs,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double val = alpha*(std::pow(double(sD.time),-std::exp(power))
		      - std::pow(double(prevTime),-std::exp(power)));
  prevTime = (unsigned int)(sD.time);
  std::for_each(probs.begin(),probs.end(),
		[&val](double & x){
		  x += val;
		});
}
							   

std::vector<double> ParamTrendPow::partial(const int notNode,
					   const int infNode,
					   const SimData & sD,
					   const TrtData & tD,
					   const FixedData & fD,
					   const DynamicData & dD){
  return std::vector<double>(parsSize,sD.time);
}
