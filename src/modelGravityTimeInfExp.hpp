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


class GravityTimeInfExpModel : public BaseModel {
 public:

  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD);

  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD);
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD);

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const std::vector<double> & mPV);
  

  GravityTimeInfExpParam mP;

  GravityTimeInfExpMcmc mcmc;

  double tuneTrt(const FixedData & fD);
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
