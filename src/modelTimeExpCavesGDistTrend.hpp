#ifndef MODEL_TIME_EXP_CAVES_G_DIST_TREND_HPP__
#define MODEL_TIME_EXP_CAVES_G_DIST_TREND_HPP__


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
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpCavesTrend.hpp"


class ModelTimeExpCavesGDistTrend : public ModelBase {
 protected:
 public:
  ModelTimeExpCavesGDistTrend(){ };
  ModelTimeExpCavesGDistTrend(const FixedData & fD);
  ModelTimeExpCavesGDistTrend(const ModelTimeExpCavesGDistTrend & m);

  virtual ModelTimeExpCavesGDistTrend &
  operator=(const ModelTimeExpCavesGDistTrend & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpCavesTrendMcmc mcmc;
};


class ModelTimeExpCavesGDistTrendFitData {
 public:
  ModelTimeExpCavesGDistTrendFitData(const ModelTimeExpCavesGDistTrend & m,
				     const std::vector<double> & all,
				     const FixedData & fD,
				     const
				     std::vector<std::vector<int> > & history);

  ModelTimeExpCavesGDistTrend m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesGDistTrendFitObjFn (const gsl_vector * x, void * params);




#endif
