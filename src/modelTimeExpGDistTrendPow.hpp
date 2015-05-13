#ifndef MODEL_TIME_EXP_G_DIST_TREND_POW_HPP__
#define MODEL_TIME_EXP_G_DIST_TREND_POW_HPP__


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
#include "paramTimeExp.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpTrendPow.hpp"


class ModelTimeExpGDistTrendPow : public ModelBase {
 protected:
 public:
  ModelTimeExpGDistTrendPow(){ };
  ModelTimeExpGDistTrendPow(const FixedData & fD);
  ModelTimeExpGDistTrendPow(const ModelTimeExpGDistTrendPow & m);

  virtual ModelTimeExpGDistTrendPow &
  operator=(const ModelTimeExpGDistTrendPow & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpTrendPowMcmc mcmc;
};


class ModelTimeExpGDistTrendPowFitData {
 public:
  ModelTimeExpGDistTrendPowFitData(const ModelTimeExpGDistTrendPow & m,
				   const std::vector<double> & all,
				   const FixedData & fD,
				   const
				   std::vector<std::vector<int> > & history);

  ModelTimeExpGDistTrendPow m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpGDistTrendPowFitObjFn (const gsl_vector * x, void * params);




#endif
