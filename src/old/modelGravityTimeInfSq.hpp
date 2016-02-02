#ifndef MODEL_GRAVITY_TIME_INF_SQ_HPP__
#define MODEL_GRAVITY_TIME_INF_SQ_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfSq.hpp"
#include "mcmcGravityTimeInfSq.hpp"


class GravityTimeInfSqModel : public BaseModel<GravityTimeInfSqParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfSqParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfSqParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfSqParam & mP, const GravityTimeInfSqParam mPInit);

  GravityTimeInfSqMcmc mcmc;


  double tuneTrt(const FixedData & fD, const GravityTimeInfSqParam & gP);
};


class GravityTimeInfSqModelFitData {
 public:
  GravityTimeInfSqModelFitData(const GravityTimeInfSqModel & m,
			       const GravityTimeInfSqParam & mP,
			       const SimData & sD,
			       const FixedData & fD,
			       const std::vector<std::vector<int> > & history);

  GravityTimeInfSqModel m;
  GravityTimeInfSqParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double gravityTimeInfSqModelFitObjFn (const gsl_vector * x, void * params);




#endif
