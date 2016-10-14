#ifndef GRAVITY_2_MCMC_HPP
#define GRAVITY_2_MCMC_HPP

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"
#include "densityEst.hpp"


class Gravity2Samples{
public:
  int numSamples;
  int numBurn;
  int numCovar;

  std::vector<double> intcp,beta,betaInf,alpha,power,trtPre,trtAct;
  std::vector<double> intcpBurn,betaBurn,betaInfBurn,alphaBurn,
    powerBurn,trtPreBurn,trtActBurn;

  double intcpSet;
  std::vector<double> betaSet,betaInfSet;
  double alphaSet;
  double powerSet;
  double trtPreSet;
  double trtActSet;

  std::vector<double> ll;
  std::vector<double> llBurn;
  double llPt,pD,Dbar,DIC;

  void setMean();
  void setMode();
  void setRand();

  void setPar(const int i,const bool fromBurn = false);

  std::vector<double> getPar() const;
};


class Gravity2Mcmc{
public:

  void load(const std::vector<std::vector<int> > & history,
    const std::vector<int> & status,
    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
    const FixedData & fD);

  double priorTrtMean;

  // MCMC samples
  Gravity2Samples samples;

  // history information of simulation
  int numNodes,T;
  std::vector<int> infHist;
  std::vector<int> trtPreHist;
  std::vector<int> trtActHist;
  std::vector<double> timeInf;

  // information about each county
  std::vector<double> d;
  std::vector<double> cc;
  std::vector<double> covar;
  int numCovar;

  std::vector<double> covarBeta_cur;
  std::vector<double> covarBeta_can;
  std::vector<double> covarBetaInf_cur;
  std::vector<double> covarBetaInf_can;
  std::vector<double> alphaW_cur;
  std::vector<double> alphaW_can;

  // current iteration of the parameters
  double intcp_cur;
  std::vector<double> beta_cur;
  std::vector<double> betaInf_cur;
  double alpha_cur;
  double power_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  std::vector<double> beta_can;
  std::vector<double> betaInf_can;
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
  void sample(int const numSamples, int const numBurn,
    const bool saveBurn = false);
  void sample(int const numSamples, int const numBurn,
    const std::vector<double> & par,
    const bool saveBurn = false);
  double ll();

  inline static void updateAlphaW(std::vector<double> & alphaW,
    const double & alphaOld,
    const double & alphaNew,
    const int numNodes);
  inline static void updateAlphaW(std::vector<double> & alphaW,
    const std::vector<double> & d,
    const std::vector<double> & cc,
    const double & alpha,
    const double & powerNew,
    const int numNodes);
  inline static void updateCovarBeta(std::vector<double> & covarBeta,
    const std::vector<double> & covar,
    const std::vector<double> & beta,
    const int numNodes,
    const int numCovar);
  inline static void updateCovarBeta(std::vector<double> & covarBeta,
    const std::vector<double> & covar,
    const double & betaOld,
    const double & betaNew,
    const int covarInd,
    const int numCovar);
};



inline void Gravity2Mcmc::updateAlphaW(std::vector<double> & alphaW,
  const double & alphaOld,
  const double & alphaNew,
  const int numNodes){
  double scale = alphaNew/alphaOld;
  int i,j;
  for(i = 0; i < numNodes; ++i)
    for(j = i; j < numNodes; ++j)
      alphaW.at(i*numNodes + j) *= scale;
}


inline void Gravity2Mcmc::updateAlphaW(std::vector<double> & alphaW,
  const std::vector<double> & d,
  const std::vector<double> & cc,
  const double & alpha,
  const double & powerNew,
  const int numNodes){
  int i,j;
  for(i = 0; i < numNodes; ++i)
    for(j = i; j < numNodes; ++j)
      alphaW.at(i*numNodes + j) = alpha * d.at(i*numNodes + j)/
        std::pow(cc.at(i*numNodes + j),std::exp(powerNew));
}

inline void Gravity2Mcmc::updateCovarBeta(std::vector<double> & covarBeta,
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


inline void Gravity2Mcmc::updateCovarBeta(std::vector<double> & covarBeta,
  const std::vector<double> & covar,
  const double & betaOld,
  const double & betaNew,
  const int covarInd,
  const int numCovar){
  int i = 0;
  double diff = betaNew - betaOld;
  std::for_each(covarBeta.begin(),covarBeta.end(),
		[&covar,&numCovar,&covarInd,&i,&diff](double & x){
		  x += covar.at(i*numCovar + covarInd)*diff;
		  ++i;
		});
}





#endif
