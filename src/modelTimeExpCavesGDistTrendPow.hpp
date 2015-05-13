#ifndef MODEL_TIME_EXP_CAVES_G_DIST_TREND_POW_HPP__
#define MODEL_TIME_EXP_CAVES_G_DIST_TREND_POW_HPP__


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
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpCavesTrendPow.hpp"


class ModelTimeExpCavesGDistTrendPow : public ModelBase {
 protected:
 public:
  ModelTimeExpCavesGDistTrendPow(){ };
  ModelTimeExpCavesGDistTrendPow(const FixedData & fD);
  ModelTimeExpCavesGDistTrendPow(const ModelTimeExpCavesGDistTrendPow & m);

  virtual ModelTimeExpCavesGDistTrendPow &
  operator=(const ModelTimeExpCavesGDistTrendPow & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpCavesTrendPowMcmc mcmc;
};


class ModelTimeExpCavesGDistTrendPowFitData {
 public:
  ModelTimeExpCavesGDistTrendPowFitData
  (const ModelTimeExpCavesGDistTrendPow & m,
   const std::vector<double> & all,
   const FixedData & fD,
   const
   std::vector<std::vector<int> > & history);

  ModelTimeExpCavesGDistTrendPow m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesGDistTrendPowFitObjFn (const gsl_vector * x, void * params);




#endif
