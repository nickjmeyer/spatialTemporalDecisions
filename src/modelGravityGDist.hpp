#ifndef MODEL_GRAVITY_G_DIST_HPP__
#define MODEL_GRAVITY_G_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrt.hpp"


class ModelGravityGDist : public ModelBase {
 protected:
 public:
  ModelGravityGDist(){ };
  ModelGravityGDist(const FixedData & fD);
  ModelGravityGDist(const ModelGravityGDist & m);

  virtual void read();

  virtual void save() const;
  
  virtual ModelGravityGDist & operator=(const ModelGravityGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);

  double tuneTrt(const FixedData & fD);
};


class ModelGravityGDistFitData {
 public:
  ModelGravityGDistFitData(const ModelGravityGDist & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const std::vector<std::vector<int> > & history);

  ModelGravityGDist m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGravityGDistFitObjFn (const gsl_vector * x, void * params);




#endif
