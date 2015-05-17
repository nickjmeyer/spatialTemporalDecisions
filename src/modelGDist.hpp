#ifndef MODEL_G_DIST_HPP__
#define MODEL_G_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramGDist.hpp"
#include "paramTrt.hpp"
#include "mcmcGDist.hpp"


class ModelGDist : public ModelBase {
 protected:
 public:
  ModelGDist(){ };
  ModelGDist(const FixedData & fD);
  ModelGDist(const ModelGDist & m);

  virtual ModelGDist & operator=(const ModelGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  GDistMcmc mcmc;
};


class ModelGDistFitData {
 public:
  ModelGDistFitData(const ModelGDist & m,
		    const std::vector<double> & all,
		    const FixedData & fD,
		    const std::vector<std::vector<int> > & history);
  
  ModelGDist m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGDistFitObjFn (const gsl_vector * x, void * params);


#endif
