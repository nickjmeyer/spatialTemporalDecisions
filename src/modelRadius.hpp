#ifndef MODEL_RADIUS_HPP__
#define MODEL_RADIUS_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramRadius.hpp"
#include "paramTrt.hpp"
#include "mcmcRadius.hpp"


class ModelRadius : public ModelBase {
 protected:
 public:
  ModelRadius(){ };
  ModelRadius(const FixedData & fD);
  ModelRadius(const ModelRadius & m);

  virtual ModelRadius & operator=(const ModelRadius & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
  
  RadiusMcmc mcmc;
};


class ModelRadiusFitData {
 public:
  ModelRadiusFitData(const ModelRadius & m,
		     const std::vector<double> & all,
		     const FixedData & fD,
		     const std::vector<std::vector<int> > & history);
  
  ModelRadius m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelRadiusFitObjFn (const gsl_vector * x, void * params);


#endif
