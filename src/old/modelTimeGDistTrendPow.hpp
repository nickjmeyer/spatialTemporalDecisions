#ifndef MODEL_TIME_G_DIST_TREND_POW_HPP__
#define MODEL_TIME_G_DIST_TREND_POW_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "timer.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrendPow.hpp"
#include "paramTime.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfTrendPow.hpp"


class ModelTimeGDistTrendPow : public ModelBase {
 protected:
 public:
  ModelTimeGDistTrendPow(){ };
  ModelTimeGDistTrendPow(const FixedData & fD);
  ModelTimeGDistTrendPow(const ModelTimeGDistTrendPow & m);

  virtual ModelTimeGDistTrendPow & operator=(const ModelTimeGDistTrendPow & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfTrendPowMcmc mcmc;
};


class ModelTimeGDistTrendPowFitData {
 public:
  ModelTimeGDistTrendPowFitData(const ModelTimeGDistTrendPow & m,
				const std::vector<double> & all,
				const FixedData & fD,
				const
				std::vector<std::vector<int> > & history);

  ModelTimeGDistTrendPow m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeGDistTrendPowFitObjFn (const gsl_vector * x, void * params);




#endif
