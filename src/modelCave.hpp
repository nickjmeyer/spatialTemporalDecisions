#ifndef MODEL_CAVES_HPP__
#define MODEL_CAVES_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "modelParamCave.hpp"
#include "mcmcCave.hpp"


class CaveModel {
 public:
  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD,
		    CaveParam & mP) const;

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const CaveParam & mP) const;
  
  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			CaveParam & mP) const;
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      CaveParam & mP);

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   CaveParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   CaveParam & mP, const CaveParam & mPInit);

  CaveMcmc mcmc;

  Estimation fitType;
};


class CaveModelFitData {
 public:
  CaveModelFitData(const CaveModel & m, const CaveParam & mP,
		      const SimData & sD,
		      const FixedData & fD,
		      const std::vector<std::vector<int> > & history);

  CaveModel m;
  CaveParam mP;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double CaveModelFitObjFn (const gsl_vector * x, void * params);


#endif
