#ifndef MODEL_TIME_EXP_CAVES_HPP__
#define MODEL_TIME_EXP_CAVES_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "timer.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravity.hpp"
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpCaves.hpp"


class ModelTimeExpCaves : public ModelBase {
 protected:
 public:
  ModelTimeExpCaves(){ };
  ModelTimeExpCaves(const FixedData & fD);
  ModelTimeExpCaves(const ModelTimeExpCaves & m);

  virtual void read();

  virtual ModelTimeExpCaves & operator=(const ModelTimeExpCaves & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);

  double tuneTrt(const FixedData & fD);
  
  GravityTimeInfExpCavesMcmc mcmc;
};


class ModelTimeExpCavesFitData {
 public:
  ModelTimeExpCavesFitData(const ModelTimeExpCaves & m,
			   const std::vector<double> & all,
			   const FixedData & fD,
			   const
			   std::vector<std::vector<int> > & history);

  ModelTimeExpCaves m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesFitObjFn (const gsl_vector * x, void * params);




#endif
