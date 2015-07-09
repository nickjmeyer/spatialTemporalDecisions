#ifndef GDIST_MCMC_HPP__
#define GDIST_MCMC_HPP__

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"


class GDistSamples{
 public:
  int numSamples;
  int numBurn;
  int numCovar;

  std::vector<double> intcp,beta,alpha,trtPre,trtAct;
  std::vector<double> intcpHist,betaHist,alphaHist,trtPreHist,trtActHist;

  double intcpSet;
  std::vector<double> betaSet;
  double alphaSet;
  double trtPreSet;
  double trtActSet;

  std::vector<double> ll;
  std::vector<double> llHist;
  double llPt,pD,Dbar,DIC;

  void setMean();
  void setRand();

  void setPar(const int i,const bool fromBurn = false);

  std::vector<double> getPar() const;
};


class GDistMcmc{
 public:

  void load(const std::vector<std::vector<int> > & history,
	    const std::vector<int> & status,
	    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
	    const FixedData & fD);

  double priorTrtMean;

  // MCMC samples
  GDistSamples samples;

  // history information of simulation
  int numNodes,T;
  std::vector<int> infHist;
  std::vector<int> trtPreHist;
  std::vector<int> trtActHist;
  std::vector<double> timeInf;

  // information about each county
  std::vector<double> d;
  std::vector<double> covar;
  int numCovar;

  // current iteration of the parameters
  double intcp_cur;
  double alpha_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  double alpha_can;
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
  void sample(int const numSamples, int const numBurn,
	      const bool saveBurn = false);
  void sample(int const numSamples, int const numBurn,
	      const std::vector<double> & par,
	      const bool saveBurn = false);
  double ll();

};




#endif
