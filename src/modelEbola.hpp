#ifndef MODEL_EBOLA_HPP__
#define MODEL_EBOLA_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "modelParamEbola.hpp"


class EbolaModel {
 public:
  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD,
		    EbolaParam & mP) const;

  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const EbolaParam & mP) const;
  
  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			EbolaParam & mP) const;
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      EbolaParam & mP);

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   EbolaParam & mP);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   EbolaParam & mP, const EbolaParam & mPInit);
};


class EbolaModelFitData {
 public:
  EbolaModelFitData(const EbolaModel & m, const EbolaParam & mP,
		      const SimData & sD,
		      const FixedData & fD,
		      const std::vector<std::vector<int> > & history);

  EbolaModel m;
  EbolaParam mP;
  FixedData fD;
  std::vector< std::vector<int> > history;
};


double ebolaModelFitObjFn (const gsl_vector * x, void * params);


#endif
