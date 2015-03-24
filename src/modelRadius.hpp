#ifndef MODEL_RADIUS_HPP__
#define MODEL_RADIUS_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "modelParamRadius.hpp"
#include "mcmcRadius.hpp"


class RadiusModel {
 public:
  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD,
		    RadiusParam & mP) const;

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const RadiusParam & mP) const;
  
  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			RadiusParam & mP) const;
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      RadiusParam & mP);

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   RadiusParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   RadiusParam & mP, const RadiusParam mPInit);

  RadiusMcmc mcmc;

  Estimation fitType;
};


class RadiusModelFitData {
 public:
  RadiusModelFitData(const RadiusModel & m, const RadiusParam & mP,
		     const SimData & sD,
		     const FixedData & fD,
		     const std::vector<std::vector<int> > & history);

  RadiusModel m;
  RadiusParam mP;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double RadiusModelFitObjFn (const gsl_vector * x, void * params);


#endif
