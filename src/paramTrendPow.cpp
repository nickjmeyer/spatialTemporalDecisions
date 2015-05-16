#include "paramTrendPow.hpp"


unsigned int ParamTrendPow::initParsSize(const FixedData & fD){
  return 2U; // unsigned int literal
}


std::vector<std::string> ParamTrendPow::initNames(){
  return {"trend","trendPow"};
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
  double val = alpha*std::pow(double(sD.time+1),power);
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
  double val = alpha*(std::pow(double(sD.time+1),power)
		      - std::pow(double(prevTime+1),power));
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
  double alpha = pars.at(0);
  double power = pars.at(1);
  double tPower = std::pow(double(sD.time+1),power);
  
  return {tPower,alpha*std::log(double(sD.time+1))*tPower};
}



std::vector<double> ParamTrendPow::partial2(const int notNode,
					    const int infNode,
					    const SimData & sD,
					    const TrtData & tD,
					    const FixedData & fD,
					    const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double tPower = std::pow(double(sD.time+1),power);
  double logT = std::log(double(sD.time+1));
  
  return {0.0, logT*tPower, logT*tPower, alpha*logT*logT*tPower};
}
