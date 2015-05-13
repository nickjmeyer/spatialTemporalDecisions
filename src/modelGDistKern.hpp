#ifndef MODEL_G_DIST_KERN_HPP__
#define MODEL_G_DIST_KERN_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramGDistKern.hpp"
#include "paramTrt.hpp"


class ModelGDistKern : public ModelBase {
 protected:
 public:
  ModelGDistKern() { };
  ModelGDistKern(const FixedData & fD);
  ModelGDistKern(const ModelGDistKern & m);

  virtual ModelBase * clone() const {return new ModelGDistKern(*this);};

  virtual ModelGDistKern & operator=(const ModelGDistKern & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class ModelGDistKernFitData {
 public:
  ModelGDistKernFitData(const ModelGDistKern & m,
			const std::vector<double> & all,
			const FixedData & fD,
			const std::vector<std::vector<int> > & history);
  
  ModelGDistKern m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGDistKernFitObjFn (const gsl_vector * x, void * params);


#endif
