#ifndef MODEL_HPP__
#define MODEL_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "param.hpp"

enum Estimation {INVALID = -1, MLE = 0, MLES = 1, MCMC = 2};

class ModelBase {
 protected:
  unsigned int numPars;
  int set;
  std::vector<ParamBase *> pars;
  std::vector<double> probs;
  std::vector<double> expitInfProbs;
  std::vector<double> expitRevProbs;
  std::vector<double> quick;
  std::vector<double> fisher;
  Eigen::VectorXd meanHit;
  Eigen::MatrixXd varHit;
  int ready;
  int numInfected,numNotInfec;

 public:
  ModelBase(){ };
  ModelBase(const std::vector<ParamBase *> & newPars,
	    const FixedData & fD);
  virtual ~ModelBase();

  virtual void read() {
    std::cout << "Error: ModelBase::read(): model is not setup to save"
	      << std::endl;
    throw(1);
  }


  virtual void save() const {
    std::cout << "Error: ModelBase::save(): model is not setup to save"
	      << std::endl;
    throw(1);
  }

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

  virtual void setQuick(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD);

  virtual std::vector<double> & getQuick();
  
  virtual double oneOnOne(const int notNode, const int infNode,
			  const int numNodes) const;

  virtual std::vector<double> partial(const int notNode,
				      const int infNode,
				      const SimData & sD,
				      const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD);

  virtual std::vector<double> partial2(const int notNode,
				       const int infNode,
				       const SimData & sD,
				       const TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD);

  virtual void setFisher(const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD);

  virtual bool sample();
  
  virtual void revert();

  virtual std::vector<double> getFisher() {return fisher;};

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
