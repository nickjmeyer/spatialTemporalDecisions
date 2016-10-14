#ifndef PARAM_BETA_HPP
#define PARAM_BETA_HPP

#include "param.hpp"


class ParamBeta : public ParamBase {
 protected:
  std::vector<double> covarBeta;
  std::vector<double> covar;

  int numNodes;

  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual void initInternal(const FixedData & fD);

  virtual void updateBefore();

  virtual void updateAfter();

 public:
  ParamBeta() { };
  virtual ParamBase * clone() const { return new ParamBeta(*this); };

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
