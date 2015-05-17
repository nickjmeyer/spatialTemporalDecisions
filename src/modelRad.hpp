#ifndef MODEL_RAD_HPP__
#define MODEL_RAD_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGDist.hpp"
#include "paramRad.hpp"
#include "paramTrt.hpp"


class ModelRad : public ModelBase {
 protected:
 public:
  ModelRad(){ };
  ModelRad(const FixedData & fD);
  ModelRad(const ModelRad & m);

  virtual ModelRad & operator=(const ModelRad & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class ModelRadFitData {
 public:
  ModelRadFitData(const ModelRad & m,
		  const std::vector<double> & all,
		  const FixedData & fD,
		  const std::vector<std::vector<int> > & history);
  
  ModelRad m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelRadFitObjFn (const gsl_vector * x, void * params);


#endif
