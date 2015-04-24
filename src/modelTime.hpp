#ifndef MODEL_TIME_HPP__
#define MODEL_TIME_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "timer.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravity.hpp"
#include "paramTime.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInf.hpp"


class ModelTime : public ModelBase {
 protected:
 public:
  ModelTime(){ };
  ModelTime(const FixedData & fD);
  ModelTime(const ModelTime & m);

  virtual void read();

  virtual ModelTime & operator=(const ModelTime & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);
  
  GravityTimeInfMcmc mcmc;
};


class ModelTimeFitData {
 public:
  ModelTimeFitData(const ModelTime & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const
			   std::vector<std::vector<int> > & history);

  ModelTime m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeFitObjFn (const gsl_vector * x, void * params);




#endif
