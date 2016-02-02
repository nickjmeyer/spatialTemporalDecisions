#ifndef MODEL_TIME_G_DIST_TREND_HPP__
#define MODEL_TIME_G_DIST_TREND_HPP__


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
#include "paramTrend.hpp"
#include "paramTime.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfTrend.hpp"


class ModelTimeGDistTrend : public ModelBase {
 protected:
 public:
  ModelTimeGDistTrend(){ };
  ModelTimeGDistTrend(const FixedData & fD);
  ModelTimeGDistTrend(const ModelTimeGDistTrend & m);

  virtual ModelTimeGDistTrend & operator=(const ModelTimeGDistTrend & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfTrendMcmc mcmc;
};


class ModelTimeGDistTrendFitData {
 public:
  ModelTimeGDistTrendFitData(const ModelTimeGDistTrend & m,
			     const std::vector<double> & all,
			     const FixedData & fD,
			     const
			     std::vector<std::vector<int> > & history);

  ModelTimeGDistTrend m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeGDistTrendFitObjFn (const gsl_vector * x, void * params);




#endif
