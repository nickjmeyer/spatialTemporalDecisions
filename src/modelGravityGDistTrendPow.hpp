#ifndef MODEL_GRAVITY_G_DIST_TREND_POW_HPP__
#define MODEL_GRAVITY_G_DIST_TREND_POW_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrendPow.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTrendPow.hpp"


class ModelGravityGDistTrendPow : public ModelBase {
 protected:
 public:
  ModelGravityGDistTrendPow(){ };
  ModelGravityGDistTrendPow(const FixedData & fD);
  ModelGravityGDistTrendPow(const ModelGravityGDistTrendPow & m);

  virtual ModelGravityGDistTrendPow &
  operator=(const ModelGravityGDistTrendPow & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  double tuneTrt(const FixedData & fD);

  GravityTrendPowMcmc mcmc;
};


class ModelGravityGDistTrendPowFitData {
 public:
  ModelGravityGDistTrendPowFitData
  (const ModelGravityGDistTrendPow & m,
   const std::vector<double> & all,
   const FixedData & fD,
   const std::vector<std::vector<int> > & history);

  ModelGravityGDistTrendPow m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravityGDistTrendPowFitObjFn (const gsl_vector * x, void * params);




#endif
