#ifndef MODEL_DIST_KERN_HPP__
#define MODEL_DIST_KERN_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramDistKern.hpp"
#include "paramTrt.hpp"


class ModelDistKern : public ModelBase {
 protected:
 public:
  ModelDistKern();
  ModelDistKern(const FixedData & fD);
  ModelDistKern(const ModelDistKern & m);

  virtual ModelBase * clone() const {return new ModelDistKern(*this);};

  virtual void read();
  
  virtual ModelDistKern & operator=(const ModelDistKern & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class ModelDistKernFitData {
 public:
  ModelDistKernFitData(const ModelDistKern & m,
		       const std::vector<double> & all,
		       const FixedData & fD,
		       const std::vector<std::vector<int> > & history);
  
  ModelDistKern m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelDistKernFitObjFn (const gsl_vector * x, void * params);


#endif
