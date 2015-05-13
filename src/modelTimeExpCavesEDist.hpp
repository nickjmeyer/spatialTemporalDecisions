#ifndef MODEL_TIME_EXP_CAVES_E_DIST_HPP__
#define MODEL_TIME_EXP_CAVES_E_DIST_HPP__


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
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"

class ModelTimeExpCavesEDist : public ModelBase {
 protected:
 public:
  ModelTimeExpCavesEDist(){ };
  ModelTimeExpCavesEDist(const FixedData & fD);
  ModelTimeExpCavesEDist(const ModelTimeExpCavesEDist & m);

  virtual ModelTimeExpCavesEDist & operator=(const ModelTimeExpCavesEDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);
};


class ModelTimeExpCavesEDistFitData {
 public:
  ModelTimeExpCavesEDistFitData(const ModelTimeExpCavesEDist & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const
			   std::vector<std::vector<int> > & history);

  ModelTimeExpCavesEDist m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesEDistFitObjFn (const gsl_vector * x, void * params);




#endif
