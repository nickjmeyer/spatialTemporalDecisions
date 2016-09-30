#ifndef EDGE_TO_EDGE_2_MCMC_HPP__
#define EDGE_TO_EDGE_2_MCMC_HPP__

#include <vector>
#include <algorithm>
#include "rand.hpp"
#include "data.hpp"
#include "densityEst.hpp"


class EdgeToEdge2Samples{
public:
  int numSamples;
  int numBurn;
  int numCovar;

  std::vector<double> intcp,beta,betaInf,trtPre,trtAct;
  std::vector<double> intcpBurn,betaBurn,betaInfBurn,
    trtPreBurn,trtActBurn;

  double intcpSet;
  std::vector<double> betaSet,betaInfSet;
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


class EdgeToEdge2Mcmc{
public:

  void load(const std::vector<std::vector<int> > & history,
    const std::vector<int> & status,
    const FixedData & fD);
  void load(const std::vector<std::vector<int> > & history,
    const FixedData & fD);

  double priorTrtMean;

  // MCMC samples
  EdgeToEdge2Samples samples;

  // history information of simulation
  int numNodes,T;
  std::vector<int> infHist;
  std::vector<int> trtPreHist;
  std::vector<int> trtActHist;
  std::vector<double> timeInf;

  // information about each county
  std::vector<int> net;
  std::vector<double> d;
  std::vector<double> cc;
  std::vector<double> covar;
  int numCovar;

  std::vector<double> covarBeta_cur;
  std::vector<double> covarBeta_can;
  std::vector<double> covarBetaInf_cur;
  std::vector<double> covarBetaInf_can;

  // current iteration of the parameters
  double intcp_cur;
  std::vector<double> beta_cur;
  std::vector<double> betaInf_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  std::vector<double> beta_can;
  std::vector<double> betaInf_can;
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



inline void EdgeToEdge2Mcmc::updateCovarBeta(std::vector<double> & covarBeta,
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


inline void EdgeToEdge2Mcmc::updateCovarBeta(std::vector<double> & covarBeta,
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
