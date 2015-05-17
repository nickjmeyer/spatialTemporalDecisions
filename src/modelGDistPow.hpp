#ifndef MODEL_G_DIST_POW_HPP__
#define MODEL_G_DIST_POW_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramGDistPow.hpp"
#include "paramTrt.hpp"
#include "mcmcGDistPow.hpp"


class ModelGDistPow : public ModelBase {
 protected:
 public:
  ModelGDistPow(){ };
  ModelGDistPow(const FixedData & fD);
  ModelGDistPow(const ModelGDistPow & m);

  virtual ModelGDistPow & operator=(const ModelGDistPow & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  GDistPowMcmc mcmc;
};


class ModelGDistPowFitData {
 public:
  ModelGDistPowFitData(const ModelGDistPow & m,
		       const std::vector<double> & all,
		       const FixedData & fD,
		       const std::vector<std::vector<int> > & history);
  
  ModelGDistPow m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGDistPowFitObjFn (const gsl_vector * x, void * params);


#endif
