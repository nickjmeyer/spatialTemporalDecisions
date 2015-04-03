#ifndef MODEL_GRAVITY_TIME_INF_HPP__
#define MODEL_GRAVITY_TIME_INF_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamGravityTimeInf.hpp"
#include "mcmcGravityTimeInf.hpp"


class GravityTimeInfModel : public BaseModel {
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

  GravityTimeInfParam mP;
  
  GravityTimeInfMcmc mcmc;


  double tuneTrt(const FixedData & fD);
};


class GravityTimeInfModelFitData {
 public:
  GravityTimeInfModelFitData(const GravityTimeInfModel & m,
			     const GravityTimeInfParam & mP,
			     const SimData & sD,
			     const FixedData & fD,
			     const std::vector<std::vector<int> > & history);

  GravityTimeInfModel m;
  GravityTimeInfParam mP;
  SimData sD;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double gravityTimeInfModelFitObjFn (const gsl_vector * x, void * params);




#endif
