#ifndef MODEL_GRAVITY_HPP__
#define MODEL_GRAVITY_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravity.hpp"
#include "paramTrt.hpp"
#include "mcmc.hpp"


class ModelGravity : public ModelBase {
 protected:
 public:
  ModelGravity(){ };
  ModelGravity(const FixedData & fD);
  ModelGravity(const ModelGravity & m);

  virtual void read();
  
  virtual ModelGravity & operator=(const ModelGravity & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  double tuneTrt(const FixedData & fD);
  
  GravityMcmc mcmc;
};


class ModelGravityFitData {
 public:
  ModelGravityFitData(const ModelGravity & m,
		      const std::vector<double> & all,
		      const SimData & sD,
		      const FixedData & fD,
		      const std::vector<std::vector<int> > & history);

  ModelGravity m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravityFitObjFn (const gsl_vector * x, void * params);




#endif
