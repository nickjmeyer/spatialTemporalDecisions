#include "osspAgent.hpp"

template class OsspAgent<GravityTimeInfExpCavesModel>;


OsspAgentTuneParam::OsspAgentTuneParam(){
  eps = 0.8;
}



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
  int ind = 0;
  if(qvalues.size() > 1){
    int i,I;
    double mx = *std::max_element(qvalues.begin(),qvalues.end());
    
    if(njm::runif01() < tp.eps){
      // uniformly take best actions
      std::vector<int> mInd;
      I = qvalues.size();
      for(i = 0; i < I; ++i){
	if(qvalues.at(i) == mx){
	  mInd.push_back(i);
	}
      }

      std::priority_queue<std::pair<double,int> > pq;
      I = mInd.size();
      for(i = 0; i < I; ++i)
	pq.push(std::pair<double,int>(njm::runif01(),i));
      ind = mInd.at(pq.top().second);
    }
    else{
      // take soft max of sub best actions
      std::vector<double> mQ;
      std::vector<int> mInd;
      I = qvalues.size();
      for(i = 0; i < I; ++i){
	if(qvalues.at(i) < mx){
	  mQ.push_back(qvalues.at(i));
	  mInd.push_back(i);
	}
      }

      std::vector<double> probs = mQ;
      std::for_each(probs.begin(),probs.end(),
		    [this](double & x)
		    {
		      x = std::exp(x);
		    });

      double total = std::accumulate(probs.begin(),probs.end(),0.0);
      std::for_each(probs.begin(),probs.end(),
		    [&total](double & x)
		    {
		      x /= total;
		    });
  
      double cur = probs.at(0);
      double num = njm::runif01();
      i = 0;
      while(cur < num)
	cur += probs.at(++i);
      try{
	ind = mInd.at(i);
      }
      catch(...){
	std::cout << "i: " << i
		  << std::endl
		  << "mInd: " << mInd.size()
		  << std::endl
		  << "probs: " << probs.size()
		  << std::endl;
      }
	  
    }
  }
  else{
    ind = 0;
  }
  
  tD.a = aCand.at(ind);
  tD.p = pCand.at(ind);
}
