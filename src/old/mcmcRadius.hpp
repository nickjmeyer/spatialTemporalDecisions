#ifndef RADIUS_MCMC_HPP__
#define RADIUS_MCMC_HPP__

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"
#include "densityEst.hpp"


class RadiusSamples{
 public:
  int numSamples;
int numBurn;
  int numCovar;
  
  std::vector<double> intcp,radius,trtPre,trtAct;
std::vector<double> intcpBurn,radiusBurn,trtPreBurn,trtActBurn;

  double intcpSet;
  double radiusSet;
  double trtPreSet;
  double trtActSet;

  std::vector<double> ll;
std::vector<double> llBurn;
  double llPt,pD,Dbar,DIC;

  void setMean();
  void setMode();
  void setRand();

  std::vector<double> getPar() const;
};


class RadiusMcmc{
 public:

  void load(const std::vector<std::vector<int> > & history,
	    const std::vector<int> & status,
	    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
	    const FixedData & fD);

  double priorTrtMean;
  
  // MCMC samples
  RadiusSamples samples;

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
  double radius_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  double radius_can;
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
