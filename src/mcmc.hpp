#ifndef MCMC_HPP__
#define MCMC_HPP__

#include <armadillo>
#include "rand.hpp"
#include "data.hpp"
// #include "model.hpp"


class GravitySamples{
 public:
  arma::mat beta;
  arma::colvec intcp,alpha,power,trtPre,trtAct;

  double intcpMean;
  arma::colvec betaMean;
  double alphaMean;
  double powerMean;
  double trtPreMean;
  double trtActMean;

  arma::colvec ll;
  double llPt,pD,Dbar,DIC;
};


class GravityMcmc{
 public:

  void load(const SimData & sD, const TrtData & tD,
	    const FixedData & fD, const DynamicData & dD);
  
  // MCMC samples
  GravitySamples samples;

  // history information of simulation
  int numNodes,T;
  arma::Mat<int> infHist;
  arma::Mat<int> trtPreHist;
  arma::Mat<int> trtActHist;
  arma::Mat<double> timeInf;

  // information about each county
  int covar;
  arma::mat d;
  arma::mat cc;
  arma::mat Xcov;
  arma::colvec XcovBeta_cur;
  arma::colvec XcovBeta_can;
  arma::mat alphaW_cur;
  arma::mat alphaW_can;

  // current iteration of the parameters
  double intcp_cur;
  arma::colvec beta_cur;
  double alpha_cur;
  double power_cur;
  double trtPre_cur;
  double trtAct_cur;

  double intcp_can;
  arma::colvec beta_can;
  double alpha_can;
  double power_can;
  double trtPre_can;
  double trtAct_can;

  // likelihood values
  double ll_cur;
  double ll_can;
  arma::colvec llVec_cur;
  arma::colvec llVec_can;

  // MH step
  arma::colvec mh;
  arma::colvec acc;
  arma::colvec att;
  arma::colvec tau;

  // variables used for priors and such
  arma::colvec mu;

  //functions
  void sample(int const numSamples, int const numBurn);
  double ll(int const b, int const B);
  double ll();
  void update_alphaW(std::string const par, int const b, int const B);
  int block2Time(int const b, int const B);

};


#endif
