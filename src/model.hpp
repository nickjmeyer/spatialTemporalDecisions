#ifndef MODEL_HPP__
#define MODEL_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "modelParam.hpp"
#include "modelParamGravityTimeInf.hpp"
#include "modelParamGravityTimeInfSq.hpp"
#include "modelParamGravityTimeInfSqrt.hpp"
#include "modelParamGravityTimeInfLog.hpp"
#include "modelParamGravityTimeInfExp.hpp"
#include "modelParamGravityTimeInfExpCaves.hpp"
#include "modelParamGravityTimeInfExpLCaves.hpp"
#include "modelParamGravityTimeInfExpRCaves.hpp"
#include "mcmc.hpp"

enum Estimation {MLE = 0,MCMC = 1};

template<class MP>
class BaseModel {
 public:

  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD,
		    MP & mP) const;

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const MP & mP) const = 0;

  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			MP & mP) const;
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      MP & mP);
};



class GravityModel : public BaseModel<GravityParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityParam & mP) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   GravityParam & mP, const GravityParam mPInit);

  GravityMcmc mcmc;

  Estimation fitType;
  

  double tuneTrt(const FixedData & fD, const GravityParam & gP);
};


class GravityModelFitData {
 public:
  GravityModelFitData(const GravityModel & m, const GravityParam & mP,
		      const SimData & sD,
		      const FixedData & fD,
		      const std::vector<std::vector<int> > & history);

  GravityModel m;
  GravityParam mP;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double gravityModelFitObjFn (const gsl_vector * x, void * params);




#endif
