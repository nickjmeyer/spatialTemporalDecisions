#ifndef MODEL_TIME_EXP_CAVES_G_DIST_HPP__
#define MODEL_TIME_EXP_CAVES_G_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "timer.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"

class ModelTimeExpCavesGDist : public ModelBase {
 protected:
 public:
  ModelTimeExpCavesGDist(){ };
  ModelTimeExpCavesGDist(const FixedData & fD);
  ModelTimeExpCavesGDist(const ModelTimeExpCavesGDist & m);

  virtual void read();

  virtual void save() const;

  virtual ModelTimeExpCavesGDist & operator=(const ModelTimeExpCavesGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);
};


class ModelTimeExpCavesGDistFitData {
 public:
  ModelTimeExpCavesGDistFitData(const ModelTimeExpCavesGDist & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const
			   std::vector<std::vector<int> > & history);

  ModelTimeExpCavesGDist m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesGDistFitObjFn (const gsl_vector * x, void * params);




#endif
