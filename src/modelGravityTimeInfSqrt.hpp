#ifndef MODEL_GRAVITY_TIME_INF_SQRT_HPP__
#define MODEL_GRAVITY_TIME_INF_SQRT_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfSqrt.hpp"
#include "mcmcGravityTimeInfSqrt.hpp"


class GravityTimeInfSqrtModel : public BaseModel<GravityTimeInfSqrtParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfSqrtParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfSqrtParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfSqrtParam & mP, const GravityTimeInfSqrtParam mPInit);

  GravityTimeInfSqrtMcmc mcmc;


  double tuneTrt(const FixedData & fD, const GravityTimeInfSqrtParam & gP);
};


class GravityTimeInfSqrtModelFitData {
 public:
  GravityTimeInfSqrtModelFitData(const GravityTimeInfSqrtModel & m,
				 const GravityTimeInfSqrtParam & mP,
				 const SimData & sD,
				 const FixedData & fD,
				 const std::vector<std::vector<int> > & history);

  GravityTimeInfSqrtModel m;
  GravityTimeInfSqrtParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double gravityTimeInfSqrtModelFitObjFn (const gsl_vector * x, void * params);




#endif
