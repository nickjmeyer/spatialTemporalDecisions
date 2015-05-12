#ifndef MODEL_TIME_EXP_G_DIST_HPP__
#define MODEL_TIME_EXP_G_DIST_HPP__


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
#include "paramTimeExp.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExp.hpp"

class ModelTimeExpGDist : public ModelBase {
 protected:
 public:
  ModelTimeExpGDist(){ };
  ModelTimeExpGDist(const FixedData & fD);
  ModelTimeExpGDist(const ModelTimeExpGDist & m);

  virtual void read();

  virtual void save() const;

  virtual ModelTimeExpGDist & operator=(const ModelTimeExpGDist & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);

  GravityTimeInfExpMcmc mcmc;
};


class ModelTimeExpGDistFitData {
 public:
  ModelTimeExpGDistFitData(const ModelTimeExpGDist & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const
			   std::vector<std::vector<int> > & history);

  ModelTimeExpGDist m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpGDistFitObjFn (const gsl_vector * x, void * params);




#endif
