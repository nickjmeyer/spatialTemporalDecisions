#ifndef MODEL_GRAVITY_TIME_INF_EXP_L_CAVES_HPP__
#define MODEL_GRAVITY_TIME_INF_EXP_L_CAVES_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfExpLCaves.hpp"
#include "mcmcGravityTimeInfExpLCaves.hpp"


class
GravityTimeInfExpLCavesModel : public BaseModel<GravityTimeInfExpLCavesParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfExpLCavesParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpLCavesParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpLCavesParam & mP,
	   const GravityTimeInfExpLCavesParam mPInit);

  GravityTimeInfExpLCavesMcmc mcmc;

  Estimation fitType;
  

  double tuneTrt(const FixedData & fD, const GravityTimeInfExpLCavesParam & gP);
};


class GravityTimeInfExpLCavesModelFitData {
 public:
  GravityTimeInfExpLCavesModelFitData(const GravityTimeInfExpLCavesModel & m,
				      const GravityTimeInfExpLCavesParam & mP,
				      const SimData & sD,
				      const FixedData & fD,
				      const
				      std::vector<std::vector<int> > & history);

  GravityTimeInfExpLCavesModel m;
  GravityTimeInfExpLCavesParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
gravityTimeInfExpLCavesModelFitObjFn (const gsl_vector * x, void * params);




#endif
