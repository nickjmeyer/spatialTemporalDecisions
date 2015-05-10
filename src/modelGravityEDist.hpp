#ifndef MODEL_GRAVITY_E_DIST_HPP__
#define MODEL_GRAVITY_E_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityEDist.hpp"
#include "paramTrt.hpp"
#include "mcmc.hpp"


class ModelGravityEDist : public ModelBase {
 protected:
 public:
  ModelGravityEDist(){ };
  ModelGravityEDist(const FixedData & fD);
  ModelGravityEDist(const ModelGravityEDist & m);

  virtual void read();

  virtual void save() const;
  
  virtual ModelGravityEDist & operator=(const ModelGravityEDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class ModelGravityEDistFitData {
 public:
  ModelGravityEDistFitData(const ModelGravityEDist & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const std::vector<std::vector<int> > & history);

  ModelGravityEDist m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravityEDistFitObjFn (const gsl_vector * x, void * params);




#endif
