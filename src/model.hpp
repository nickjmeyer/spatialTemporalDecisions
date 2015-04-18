#ifndef MODEL_HPP__
#define MODEL_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "param.hpp"

enum Estimation {INVALID = -1, MLE = 0, MCMC = 1};

class ModelBase {
 protected:
  int set;
  std::vector<ParamBase *> pars;
  std::vector<double> probs;
  std::vector<double> expitInfProbs;
  std::vector<double> expitRevProbs;

 public:
  ModelBase(){ };
  ModelBase(const std::vector<ParamBase *> & newPars,
	    const FixedData & fD);
  virtual ~ModelBase();

  virtual void read() = 0;

  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD);

  virtual std::vector<double> infProbs();
  
  virtual void revProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD);

  virtual std::vector<double> revProbs();
  
  virtual void setFill(const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD);

  virtual void modFill(const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD);
  
  virtual double oneOnOne(const int notNode, const int infNode,
			  const int numNodes) const;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit) = 0;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars) = 0;
  
  virtual std::vector<double> getPar() const;
  
  virtual std::vector<double>::const_iterator
  putPar(std::vector<double>::const_iterator it);

  virtual void setType(const Estimation & est);
  virtual Estimation getType() const;
  virtual Estimation & getType();

  Estimation fitType;
};




#endif
