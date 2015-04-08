#include "osspAgent.hpp"

template class OsspAgent<GravityTimeInfExpCavesModel>;


std::vector<double> OsspAgentTuneParam::getPar() const {
  return std::vector<double> (0);
}


void OsspAgentTuneParam::putPar(const std::vector<double> & par){
}


template <class M>
OsspAgent<M>::OsspAgent(){
  name = "ossp";
}


template <class M>
void OsspAgent<M>::reset(){
}



template <class M>
void OsspAgent<M>::applyTrt(const SimData & sD,
			    TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD,
			    M & m){
  std::vector<double> probs = qvalues;
  std::for_each(probs.begin(),probs.end(),
		[this](double & x)
		{
		  x = std::exp(tp.lambda*x);
		});

  double total = std::accumulate(probs.begin(),probs.end(),0.0);
  std::for_each(probs.begin(),probs.end(),
		[&total](double & x)
		{
		  x /= total;
		});

  double cur = probs.at(0);
  double num = njm::runif01();
  int ind = 0;
  while(cur < num)
    cur += probs.at(++ind);
  tD.a = aCand.at(ind);
  tD.p = pCand.at(ind);
}
