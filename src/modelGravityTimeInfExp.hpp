#ifndef MODEL_GRAVITY_TIME_INF_EXP_HPP__
#define MODEL_GRAVITY_TIME_INF_EXP_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfExp.hpp"
#include "mcmcGravityTimeInfExp.hpp"


class GravityTimeInfExpModel : public BaseModel<GravityTimeInfExpParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfExpParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpParam & mP, const GravityTimeInfExpParam mPInit);

  GravityTimeInfExpMcmc mcmc;

  Estimation fitType;
  

  double tuneTrt(const FixedData & fD, const GravityTimeInfExpParam & gP);
};


class GravityTimeInfExpModelFitData {
 public:
  GravityTimeInfExpModelFitData(const GravityTimeInfExpModel & m,
				const GravityTimeInfExpParam & mP,
				const SimData & sD,
				const FixedData & fD,
				const
				std::vector<std::vector<int> > & history);

  GravityTimeInfExpModel m;
  GravityTimeInfExpParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double gravityTimeInfExpModelFitObjFn (const gsl_vector * x, void * params);




#endif
