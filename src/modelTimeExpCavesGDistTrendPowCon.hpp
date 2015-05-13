#ifndef MODEL_TIME_EXP_CAVES_G_DIST_TREND_POW_CON_HPP__
#define MODEL_TIME_EXP_CAVES_G_DIST_TREND_POW_CON_HPP__


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
#include "paramTrendPowCon.hpp"
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpCavesTrendPowCon.hpp"


class ModelTimeExpCavesGDistTrendPowCon : public ModelBase {
 protected:
 public:
  ModelTimeExpCavesGDistTrendPowCon(){ };
  ModelTimeExpCavesGDistTrendPowCon(const FixedData & fD);
  ModelTimeExpCavesGDistTrendPowCon
  (const ModelTimeExpCavesGDistTrendPowCon & m);

  virtual ModelTimeExpCavesGDistTrendPowCon &
  operator=(const ModelTimeExpCavesGDistTrendPowCon & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpCavesTrendPowConMcmc mcmc;
};


class ModelTimeExpCavesGDistTrendPowConFitData {
 public:
  ModelTimeExpCavesGDistTrendPowConFitData
  (const ModelTimeExpCavesGDistTrendPowCon & m,
   const std::vector<double> & all,
   const FixedData & fD,
   const
   std::vector<std::vector<int> > & history);

  ModelTimeExpCavesGDistTrendPowCon m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesGDistTrendPowConFitObjFn (const gsl_vector * x, void * params);




#endif
