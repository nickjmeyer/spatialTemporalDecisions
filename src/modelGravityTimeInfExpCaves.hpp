#ifndef MODEL_GRAVITY_TIME_INF_EXP_CAVES_HPP__
#define MODEL_GRAVITY_TIME_INF_EXP_CAVES_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfExpCaves.hpp"
#include "mcmcGravityTimeInfExpCaves.hpp"


class
GravityTimeInfExpCavesModel : public BaseModel<GravityTimeInfExpCavesParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfExpCavesParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpCavesParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpCavesParam & mP,
	   const GravityTimeInfExpCavesParam mPInit);

  GravityTimeInfExpCavesMcmc mcmc;

  Estimation fitType;
  

  double tuneTrt(const FixedData & fD, const GravityTimeInfExpCavesParam & gP);
};


class GravityTimeInfExpCavesModelFitData {
 public:
  GravityTimeInfExpCavesModelFitData(const GravityTimeInfExpCavesModel & m,
				     const GravityTimeInfExpCavesParam & mP,
				     const SimData & sD,
				     const FixedData & fD,
				     const
				     std::vector<std::vector<int> > & history);

  GravityTimeInfExpCavesModel m;
  GravityTimeInfExpCavesParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
gravityTimeInfExpCavesModelFitObjFn (const gsl_vector * x, void * params);




#endif
