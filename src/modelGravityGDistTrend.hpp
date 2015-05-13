#ifndef MODEL_GRAVITY_G_DIST_TREND_HPP__
#define MODEL_GRAVITY_G_DIST_TREND_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrend.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTrend.hpp"


class ModelGravityGDistTrend : public ModelBase {
 protected:
 public:
  ModelGravityGDistTrend(){ };
  ModelGravityGDistTrend(const FixedData & fD);
  ModelGravityGDistTrend(const ModelGravityGDistTrend & m);

  virtual ModelGravityGDistTrend & operator=(const ModelGravityGDistTrend & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  double tuneTrt(const FixedData & fD);

  GravityTrendMcmc mcmc;
};


class ModelGravityGDistTrendFitData {
 public:
  ModelGravityGDistTrendFitData(const ModelGravityGDistTrend & m,
				const std::vector<double> & all,
				const FixedData & fD,
				const std::vector<std::vector<int> > & history);

  ModelGravityGDistTrend m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravityGDistTrendFitObjFn (const gsl_vector * x, void * params);




#endif
