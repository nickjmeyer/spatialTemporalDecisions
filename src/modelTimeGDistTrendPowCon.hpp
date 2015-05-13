#ifndef MODEL_TIME_G_DIST_TREND_POW_CON_HPP__
#define MODEL_TIME_G_DIST_TREND_POW_CON_HPP__


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
#include "paramTime.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfTrendPowCon.hpp"


class ModelTimeGDistTrendPowCon : public ModelBase {
 protected:
 public:
  ModelTimeGDistTrendPowCon(){ };
  ModelTimeGDistTrendPowCon(const FixedData & fD);
  ModelTimeGDistTrendPowCon(const ModelTimeGDistTrendPowCon & m);

  virtual ModelTimeGDistTrendPowCon &
  operator=(const ModelTimeGDistTrendPowCon & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfTrendPowConMcmc mcmc;
};


class ModelTimeGDistTrendPowConFitData {
 public:
  ModelTimeGDistTrendPowConFitData(const ModelTimeGDistTrendPowCon & m,
				   const std::vector<double> & all,
				   const FixedData & fD,
				   const
				   std::vector<std::vector<int> > & history);

  ModelTimeGDistTrendPowCon m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeGDistTrendPowConFitObjFn (const gsl_vector * x, void * params);




#endif
