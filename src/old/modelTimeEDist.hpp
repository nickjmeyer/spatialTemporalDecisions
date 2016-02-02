#ifndef MODEL_TIME_E_DIST_HPP__
#define MODEL_TIME_E_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "timer.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityEDist.hpp"
#include "paramTime.hpp"
#include "paramTrt.hpp"


class ModelTimeEDist : public ModelBase {
 protected:
 public:
  ModelTimeEDist(){ };
  ModelTimeEDist(const FixedData & fD);
  ModelTimeEDist(const ModelTimeEDist & m);

  virtual ModelTimeEDist & operator=(const ModelTimeEDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);
};


class ModelTimeEDistFitData {
 public:
  ModelTimeEDistFitData(const ModelTimeEDist & m,
			const std::vector<double> & all,
			const FixedData & fD,
			const
			std::vector<std::vector<int> > & history);

  ModelTimeEDist m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeEDistFitObjFn (const gsl_vector * x, void * params);




#endif
