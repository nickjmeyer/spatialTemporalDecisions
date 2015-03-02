#ifndef CAVE_MCMC_HPP__
#define CAVE_MCMC_HPP__

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"


class CaveSamples{
 public:
  int numSamples;
  int numCovar;
  
  std::vector<double> intcp,cave,trtPre,trtAct;

  double intcpSet;
  double caveSet;
  double trtPreSet;
  double trtActSet;

  std::vector<double> ll;
  double llPt,pD,Dbar,DIC;

  void setMean();
  void setRand();

  std::vector<double> getPar() const;
};


class CaveMcmc{
 public:

  void load(const std::vector<std::vector<int> > & history,
	    const std::vector<int> & status,
	    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
	    const FixedData & fD);

  double priorTrtMean;
  
  // MCMC samples
  CaveSamples samples;

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
  

  // current iteration of the parameters
  double intcp_cur;
  double cave_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  double cave_can;
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

};




#endif
