#ifndef MODEL_RANGE_HPP__
#define MODEL_RANGE_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "modelParamRange.hpp"
#include "mcmcRange.hpp"


class RangeModel {
 public:
  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD,
		    RangeParam & mP) const;

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const RangeParam & mP) const;
  
  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			RangeParam & mP) const;
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      RangeParam & mP);

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   RangeParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   RangeParam & mP, const RangeParam mPInit);

  RangeMcmc mcmc;

  Estimation fitType;
};


class RangeModelFitData {
 public:
  RangeModelFitData(const RangeModel & m, const RangeParam & mP,
		      const SimData & sD,
		      const FixedData & fD,
		      const std::vector<std::vector<int> > & history);

  RangeModel m;
  RangeParam mP;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double RangeModelFitObjFn (const gsl_vector * x, void * params);


#endif
