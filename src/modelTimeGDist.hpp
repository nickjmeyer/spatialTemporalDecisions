#ifndef MODEL_TIME_G_DIST_HPP__
#define MODEL_TIME_G_DIST_HPP__


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
#include "paramTime.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInf.hpp"


class ModelTimeGDist : public ModelBase {
 protected:
 public:
  ModelTimeGDist(){ };
  ModelTimeGDist(const FixedData & fD);
  ModelTimeGDist(const ModelTimeGDist & m);

  virtual ModelTimeGDist & operator=(const ModelTimeGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfMcmc mcmc;
};


class ModelTimeGDistFitData {
 public:
  ModelTimeGDistFitData(const ModelTimeGDist & m,
			const std::vector<double> & all,
			const FixedData & fD,
			const
			std::vector<std::vector<int> > & history);

  ModelTimeGDist m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeGDistFitObjFn (const gsl_vector * x, void * params);




#endif
