#ifndef MODEL_2_GRAVITY_G_DIST_HPP__
#define MODEL_2_GRAVITY_G_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta2.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrt.hpp"
#include "mcmcGravity2.hpp"


class Model2GravityGDist : public ModelBase {
 protected:
 public:
  Model2GravityGDist(){ };
  Model2GravityGDist(const FixedData & fD);
  Model2GravityGDist(const Model2GravityGDist & m);

  virtual Model2GravityGDist & operator=(const Model2GravityGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  Gravity2Mcmc mcmc;
};


class Model2GravityGDistFitData {
 public:
  Model2GravityGDistFitData(const Model2GravityGDist & m,
			    const std::vector<double> & all,
			    const FixedData & fD,
			    const std::vector<std::vector<int> > & history);

  Model2GravityGDist m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravity2GDistFitObjFn (const gsl_vector * x, void * params);




#endif
