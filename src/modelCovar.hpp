#ifndef MODEL_COVAR_HPP__
#define MODEL_COVAR_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramTrt.hpp"
#include "mcmc.hpp"


class ModelCovar : public ModelBase {
 protected:
 public:
  ModelCovar(){ };
  ModelCovar(const FixedData & fD);
  ModelCovar(const ModelCovar & m);

  virtual void read();

  virtual void save() const;
  
  virtual ModelCovar & operator=(const ModelCovar & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class ModelCovarFitData {
 public:
  ModelCovarFitData(const ModelCovar & m,
		    const std::vector<double> & all,
		    const FixedData & fD,
		    const std::vector<std::vector<int> > & history);

  ModelCovar m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelCovarFitObjFn (const gsl_vector * x, void * params);




#endif
