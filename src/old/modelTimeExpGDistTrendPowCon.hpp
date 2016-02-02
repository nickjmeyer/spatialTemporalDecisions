#ifndef MODEL_TIME_EXP_G_DIST_TREND_POW_CON_HPP__
#define MODEL_TIME_EXP_G_DIST_TREND_POW_CON_HPP__


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
#include "paramTimeExp.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpTrendPowCon.hpp"


class ModelTimeExpGDistTrendPowCon : public ModelBase {
 protected:
 public:
  ModelTimeExpGDistTrendPowCon(){ };
  ModelTimeExpGDistTrendPowCon(const FixedData & fD);
  ModelTimeExpGDistTrendPowCon(const ModelTimeExpGDistTrendPowCon & m);

  virtual ModelTimeExpGDistTrendPowCon &
  operator=(const ModelTimeExpGDistTrendPowCon & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpTrendPowConMcmc mcmc;
};


class ModelTimeExpGDistTrendPowConFitData {
 public:
  ModelTimeExpGDistTrendPowConFitData(const ModelTimeExpGDistTrendPowCon & m,
				      const std::vector<double> & all,
				      const FixedData & fD,
				      const
				      std::vector<std::vector<int> > & history);

  ModelTimeExpGDistTrendPowCon m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpGDistTrendPowConFitObjFn (const gsl_vector * x, void * params);




#endif
