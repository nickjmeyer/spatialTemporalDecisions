#ifndef GRAVITY_MCMC_HPP__
#define GRAVITY_MCMC_HPP__

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"


class GravitySamples{
 public:
  int numSamples;
  int numCovar;
  
  std::vector<double> intcp,beta,alpha,power,trtPre,trtAct;

  double intcpSet;
  std::vector<double> betaSet;
  double alphaSet;
  double powerSet;
  double trtPreSet;
  double trtActSet;

  std::vector<double> ll;
  double llPt,pD,Dbar,DIC;

  void setMean();
  void setRand();

  std::vector<double> getPar() const;
};


class GravityMcmc{
 public:

  void load(const std::vector<std::vector<int> > & history,
	    const std::vector<int> & status,
	    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
	    const FixedData & fD);
  
  // MCMC samples
  GravitySamples samples;

  // history information of simulation
  int numNodes,T;
  std::vector<int> infHist;
  std::vector<int> trtPreHist;
  std::vector<int> trtActHist;
  std::vector<int> timeInf;

  // information about each county
  std::vector<double> d;
  std::vector<double> cc;
  std::vector<double> covar;
  int numCovar;
  
  std::vector<double> covarBeta_cur;
  std::vector<double> covarBeta_can;
  std::vector<double> alphaW_cur;
  std::vector<double> alphaW_can;

  // current iteration of the parameters
  double intcp_cur;
  std::vector<double> beta_cur;
  double alpha_cur;
  double power_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  std::vector<double> beta_can;
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
  double ll();

  inline static void updateAlphaW(std::vector<double> & alphaW,
				  const double & alphaOld,
				  const double & alphaNew);
  inline static void updateAlphaW(std::vector<double> & alphaW,
				  const std::vector<double> & d,
				  const std::vector<double> & cc,
				  const double & alpha,
				  const double & powerNew);
  inline static void updateCovarBeta(std::vector<double> & covarBeta,
				     const std::vector<double> & covar,
				     const std::vector<double> & beta,
				     const int numNodes,
				     const int numCovar);
};



inline void GravityMcmc::updateAlphaW(std::vector<double> & alphaW,
				      const double & alphaOld,
				      const double & alphaNew){
  double scale = alphaNew/alphaOld;
  std::for_each(alphaW.begin(),alphaW.end(),
		[&scale](double & x){x *= scale;});
}

inline void GravityMcmc::updateAlphaW(std::vector<double> & alphaW,
				      const std::vector<double> & d,
				      const std::vector<double> & cc,
				      const double & alpha,
				      const double & powerNew){
  int i = 0;
  std::for_each(alphaW.begin(),alphaW.end(),
		[&](double & x){x = alpha * d.at(i)/
		    std::pow(cc.at(i),powerNew);
		  ++i;});
}

inline void GravityMcmc::updateCovarBeta(std::vector<double> & covarBeta,
					 const std::vector<double> & covar,
					 const std::vector<double> & beta,
					 const int numNodes,
					 const int numCovar){
  int i,j;
  double prod;
  for(i = 0; i < numNodes; ++i){
    prod = 0;
    for(j = 0; j < numCovar; ++j){
      prod += covar.at(i*numCovar + j) * beta.at(j);
    }
    covarBeta.at(i) = prod;
  }
}





#endif
