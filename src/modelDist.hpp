#ifndef MODEL_DIST_HPP__
#define MODEL_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramDist.hpp"
#include "paramTrt.hpp"


class ModelDist : public ModelBase {
 protected:
 public:
  ModelDist(){ };
  ModelDist(const FixedData & fD);
  ModelDist(const ModelDist & m);

  virtual void read();
  
  virtual ModelDist & operator=(const ModelDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class ModelDistFitData {
 public:
  ModelDistFitData(const ModelDist & m,
		   const std::vector<double> & all,
		   const FixedData & fD,
		   const std::vector<std::vector<int> > & history);
  
  ModelDist m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelDistFitObjFn (const gsl_vector * x, void * params);


#endif
