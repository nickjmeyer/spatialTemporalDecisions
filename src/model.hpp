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
 // protected:
 public:
  unsigned int numPars;
  std::string name;
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
  int fixSample;

 public:
  ModelBase(){ };
  ModelBase(const std::string & str,
	    const std::vector<ParamBase *> & newPars,
	    const FixedData & fD);
  virtual ~ModelBase();

  virtual void read();

  virtual void save() const;

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

  virtual bool sample(const bool force = false);

  virtual void revert();

  virtual std::vector<double> getFisher() {return fisher;};

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit) = 0;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> pars) = 0;

  virtual std::vector<double> getPar() const;

  virtual std::vector<double>
  getPar(const std::vector<std::string> & name) const;

  virtual std::vector<double>::const_iterator
  putPar(std::vector<double>::const_iterator it);

  virtual void setPar(const std::string & name, const double & val);
  virtual void setPar(const std::vector<std::string> & name,
		      const double & val);

  virtual void linScale(const double & scale);

  virtual void setType(const Estimation & est);
  virtual void setFixSample(const int & fix);
  virtual Estimation getType() const;
  virtual Estimation & getType();

  Estimation fitType;


  virtual void estimateMle(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD);

  virtual void estimateMle(const std::vector<double> startingVals,
			   const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD);

  virtual double logll(const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD);
  virtual std::vector<double> logllGrad(const SimData & sD,
					const TrtData & tD,
					const FixedData & fD,
					const DynamicData & dD);
};


class ModelBaseFitObj {
 public:

  ModelBaseFitObj(ModelBase * const mb,
		  const SimData sD,
		  const TrtData tD,
		  const FixedData fD,
		  const DynamicData dD);

  ModelBase * mb;
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;
};


double objFn(const gsl_vector * x, void * params);

void objFnGrad(const gsl_vector * x, void * params,gsl_vector * g);

void objFnBoth(const gsl_vector * x, void * params, double * f, gsl_vector * g);






#endif
