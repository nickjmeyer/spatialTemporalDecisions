#ifndef MODEL_HPP__
#define MODEL_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "modelParam.hpp"
#include "mcmc.hpp"

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
	   GravityParam & mP, const GravityParam & mPInit);
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



class GravityModelMcmc : public BaseModel<GravityParam> {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD,
			  const GravityParam & mP) const;

  GravityMcmc mcmc;

  void sample(const SimData & sD, const TrtData & tD, const FixedData & fD);
  void sample(const SimData & sD, const TrtData & tD, const FixedData & fD,
	      const GravityParam & mP);

  void assignMean(GravityParam & mP);
  void assignMean(GravityParam & mP0, GravityParam & mP1);
  void assignRand(GravityParam & mP);
  void assignRand(GravityParam & mP0, GravityParam & mP1);
};


inline double multOneMinus(double a, double b);
  

#endif
