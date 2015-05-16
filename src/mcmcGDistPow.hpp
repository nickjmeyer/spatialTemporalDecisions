#ifndef GDIST_POW_MCMC_HPP__
#define GDIST_POW_MCMC_HPP__

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"


class GravitySamples{
 public:
  int numSamples;
  int numCovar;
  
  std::vector<double> intcp,alpha,power,trtPre,trtAct;

  double intcpSet;
  double alphaSet;
  double powerSet;
  double trtPreSet;
  double trtActSet;

  std::vector<double> ll;
  double llPt,pD,Dbar,DIC;

  void setMean();
  void setRand();

  void setPar(const int i);

  std::vector<double> getPar() const;
};


class GravityMcmc{
 public:

  void load(const std::vector<std::vector<int> > & history,
	    const std::vector<int> & status,
	    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
	    const FixedData & fD);

  double priorTrtMean;
  
  // MCMC samples
  GravitySamples samples;

  // history information of simulation
  int numNodes,T;
  std::vector<int> infHist;
  std::vector<int> trtPreHist;
  std::vector<int> trtActHist;
  std::vector<double> timeInf;

  // information about each county
  std::vector<double> d;
  int numCovar;
  
  std::vector<double> w_cur;
  std::vector<double> w_can;

  // current iteration of the parameters
  double intcp_cur;
  double alpha_cur;
  double power_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  double alpha_can;
  double power_can;
  double trtPre_can;
  double trtAct_can;

  // likelihood values
  double ll_cur;
  double ll_can;

  // MH step
  std::vector<double> mh;
  std::vector<int> acc;
  std::vector<int> att;
  // std::vector<double> tau;

  // // variables used for priors and such
  // std::vector<double> mu;

  //functions
  void sample(int const numSamples, int const numBurn);
  void sample(int const numSamples, int const numBurn,
	      const std::vector<double> & par);
  double ll();

  inline static void updateW(std::vector<double> & alphaW,
			     const double & alphaOld,
			     const double & alphaNew,
			     const int numNodes);
  inline static void updateW(std::vector<double> & alphaW,
			     const std::vector<double> & d,
			     const double & alpha,
			     const double & powerNew,
			     const int numNodes);

};


inline void GravityMcmc::updateW(std::vector<double> & w,
				 const double & alphaOld,
				 const double & alphaNew,
				 const int numNodes){
  double scale = alphaNew/alphaOld;
  int i,j;
  for(i = 0; i < numNodes; ++i)
    for(j = i; j < numNodes; ++j)
      w.at(i*numNodes + j) *= scale;
}


inline void GravityMcmc::updateW(std::vector<double> & w,
				 const std::vector<double> & d,
				 const double & alpha,
				 const double & powerNew,
				 const int numNodes){
  int i,j;
  for(i = 0; i < numNodes; ++i)
    for(j = i; j < numNodes; ++j)
      w.at(i*numNodes + j) =
	alpha * std::pow(d.at(i*numNodes + j),std::exp(powerNew));
}





#endif
