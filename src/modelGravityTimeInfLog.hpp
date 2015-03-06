#ifndef MODEL_GRAVITY_TIME_INF_LOG_HPP__
#define MODEL_GRAVITY_TIME_INF_LOG_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInfLog.hpp"
#include "mcmcGravityTimeInfLog.hpp"


class GravityTimeInfLogModel : public BaseModel<GravityTimeInfLogParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityTimeInfLogParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfLogParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityTimeInfLogParam & mP, const GravityTimeInfLogParam mPInit);

  GravityTimeInfLogMcmc mcmc;

  Estimation fitType;
  

  double tuneTrt(const FixedData & fD, const GravityTimeInfLogParam & gP);
};


class GravityTimeInfLogModelFitData {
 public:
  GravityTimeInfLogModelFitData(const GravityTimeInfLogModel & m,
				const GravityTimeInfLogParam & mP,
				const SimData & sD,
				const FixedData & fD,
				const
				std::vector<std::vector<int> > & history);

  GravityTimeInfLogModel m;
  GravityTimeInfLogParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double gravityTimeInfLogModelFitObjFn (const gsl_vector * x, void * params);




#endif
