#ifndef MODEL_HPP__
#define MODEL_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "modelParam.hpp"
#include "modelParamCave.hpp"
#include "modelParamRadius.hpp"
#include "modelParamRange.hpp"
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

class BaseModel {
 public:
  virtual ~BaseModel() {}

  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD) = 0;

  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD) = 0;
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD) = 0;

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const = 0;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit) = 0;

  virtual void setType(const Estimation & est);

  Estimation fitType;
};



class GravityModel : public BaseModel {
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
			  const DynamicData & dD) const ;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const std::vector<double> & mPV);

  GravityParam mP;

  GravityMcmc mcmc;
  
  double tuneTrt(const FixedData & fD);
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
