#ifndef MODEL_TIME_EXP_G_DIST_TREND_HPP__
#define MODEL_TIME_EXP_G_DIST_TREND_HPP__


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
#include "paramTimeExp.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpTrend.hpp"


class ModelTimeExpGDistTrend : public ModelBase {
 protected:
 public:
  ModelTimeExpGDistTrend(){ };
  ModelTimeExpGDistTrend(const FixedData & fD);
  ModelTimeExpGDistTrend(const ModelTimeExpGDistTrend & m);

  virtual ModelTimeExpGDistTrend & operator=(const ModelTimeExpGDistTrend & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpTrendMcmc mcmc;
};


class ModelTimeExpGDistTrendFitData {
 public:
  ModelTimeExpGDistTrendFitData(const ModelTimeExpGDistTrend & m,
				const std::vector<double> & all,
				const FixedData & fD,
				const
				std::vector<std::vector<int> > & history);

  ModelTimeExpGDistTrend m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpGDistTrendFitObjFn (const gsl_vector * x, void * params);




#endif
