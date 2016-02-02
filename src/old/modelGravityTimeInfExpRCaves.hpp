#ifndef MODEL_GRAVITY_TIME_INF_EXP_R_CAVES_HPP__
#define MODEL_GRAVITY_TIME_INF_EXP_R_CAVES_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfExpRCaves.hpp"
#include "mcmcGravityTimeInfExpRCaves.hpp"


class
GravityTimeInfExpRCavesModel : public BaseModel<GravityTimeInfExpRCavesParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfExpRCavesParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpRCavesParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfExpRCavesParam & mP,
	   const GravityTimeInfExpRCavesParam mPInit);

  GravityTimeInfExpRCavesMcmc mcmc;


  double tuneTrt(const FixedData & fD, const GravityTimeInfExpRCavesParam & gP);
};


class GravityTimeInfExpRCavesModelFitData {
 public:
  GravityTimeInfExpRCavesModelFitData(const GravityTimeInfExpRCavesModel & m,
				      const GravityTimeInfExpRCavesParam & mP,
				      const SimData & sD,
				      const FixedData & fD,
				      const
				      std::vector<std::vector<int> > & history);

  GravityTimeInfExpRCavesModel m;
  GravityTimeInfExpRCavesParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
gravityTimeInfExpRCavesModelFitObjFn (const gsl_vector * x, void * params);




#endif
