#ifndef MODEL_2_G_POW_G_DIST_HPP__
#define MODEL_2_G_POW_G_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta2.hpp"
#include "paramGravPowGDist.hpp"
#include "paramTrt.hpp"


class Model2GPowGDist : public ModelBase {
 protected:
 public:
  Model2GPowGDist(){ };
  Model2GPowGDist(const FixedData & fD);
  Model2GPowGDist(const Model2GPowGDist & m);

  virtual Model2GPowGDist & operator=(const Model2GPowGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars);
};


class Model2GPowGDistFitData {
 public:
  Model2GPowGDistFitData(const Model2GPowGDist & m,
			    const std::vector<double> & all,
			    const FixedData & fD,
			    const std::vector<std::vector<int> > & history);

  Model2GPowGDist m;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double modelGPow2GDistFitObjFn (const gsl_vector * x, void * params);




#endif
