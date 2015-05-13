#ifndef MODEL_GRAVITY_G_DIST_TREND_POW_CON_HPP__
#define MODEL_GRAVITY_G_DIST_TREND_POW_CON_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrendPowCon.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTrendPowCon.hpp"


class ModelGravityGDistTrendPowCon : public ModelBase {
 protected:
 public:
  ModelGravityGDistTrendPowCon(){ };
  ModelGravityGDistTrendPowCon(const FixedData & fD);
  ModelGravityGDistTrendPowCon(const ModelGravityGDistTrendPowCon & m);

  virtual ModelGravityGDistTrendPowCon &
  operator=(const ModelGravityGDistTrendPowCon & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  double tuneTrt(const FixedData & fD);

  GravityTrendPowConMcmc mcmc;
};


class ModelGravityGDistTrendPowConFitData {
 public:
  ModelGravityGDistTrendPowConFitData
  (const ModelGravityGDistTrendPowCon & m,
   const std::vector<double> & all,
   const FixedData & fD,
   const std::vector<std::vector<int> > & history);

  ModelGravityGDistTrendPowCon m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravityGDistTrendPowConFitObjFn
(const gsl_vector * x, void * params);




#endif
