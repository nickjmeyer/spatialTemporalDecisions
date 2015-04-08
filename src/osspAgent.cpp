#include "osspAgent.hpp"


std::vector<double> OsspTunePar


template <class M>
OsspAgent::OsspAgent(){
  name = "ossp";
}



template <class M>
void OsspAgent::applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & m){
  std::vector<double> probs = qvalues;
  std::for_each(probs.begin(),probs.end(),
		[](double & x)
		{
		  x = 1.0 - 1.0/(1.0 + std::exp(tp.lambda*x));
		});

  double cur = probs.at(0);
  double num = njm::runif01();
  int ind = 0;
  while(cur < num)
    cur += probs.at(+=ind);
  tD.p = pCand.at(ind);
  tD.a = aCand.at(ind);
}
