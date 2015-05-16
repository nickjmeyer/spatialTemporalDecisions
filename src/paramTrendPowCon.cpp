#include "paramTrendPowCon.hpp"


unsigned int ParamTrendPowCon::initParsSize(const FixedData & fD){
  return 2U; // unsigned int literal
}


std::vector<std::string> ParamTrendPowCon::initNames(){
  return {"trend","trendPow"};
}


void ParamTrendPowCon::initInternal(const FixedData & fD){
  prevTime = 0;
}


void ParamTrendPowCon::updateBefore(){
}


void ParamTrendPowCon::updateAfter(){
}


void ParamTrendPowCon::setFill(std::vector<double> & probs,
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


void ParamTrendPowCon::modFill(std::vector<double> & probs,
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
							   

std::vector<double> ParamTrendPowCon::partial(const int notNode,
					      const int infNode,
					      const SimData & sD,
					      const TrtData & tD,
					      const FixedData & fD,
					      const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double ePower = std::exp(power);
  double tPower = std::pow(sD.time,-ePower);
  return {tPower, -alpha*ePower*std::log(double(sD.time))*tPower};
}


std::vector<double> ParamTrendPowCon::partial2(const int notNode,
					       const int infNode,
					       const SimData & sD,
					       const TrtData & tD,
					       const FixedData & fD,
					       const DynamicData & dD){
  double alpha = pars.at(0);
  double power = pars.at(1);
  double ePower = std::exp(power);
  double tPower = std::pow(sD.time,-ePower);
  double logT = std::log(double(sD.time));
  return {0.0 , - ePower*logT*tPower,
      - ePower*logT*tPower,
      alpha*std::exp(2*power)*tPower*logT*logT - alpha*ePower*tPower*logT};
}
