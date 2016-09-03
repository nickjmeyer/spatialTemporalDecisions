#ifndef PARAM_BETA_2_HPP__
#define PARAM_BETA_2_HPP__

#include "param.hpp"


class ParamBeta2 : public ParamBase {
 protected:
  std::vector<double> covarBeta,covarBetaInf;
  std::vector<double> covar;

  int numCovar;

  int numNodes;

  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual void initInternal(const FixedData & fD);

  virtual void updateBefore();

  virtual void updateAfter();

 public:
  ParamBeta2() { };
  virtual ParamBase * clone() const { return new ParamBeta2(*this); };

  virtual void setFill(std::vector<double> & probs,
		       const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD);

  virtual void modFill(std::vector<double> & probs,
		       const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD);

  virtual std::vector<double> partial(const int notNode,
				      const int infNode,
				      const SimData & sD,
				      const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD);

};



#endif
