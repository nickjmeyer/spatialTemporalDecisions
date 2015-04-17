#ifndef PARAM_BETA_HPP__
#define PARAM_BETA_HPP__

#include "param.hpp"


class ParamBeta : public ParamBase {
 protected:
  std::vector<double> parsOld;
  std::vector<double> covar;
  std::vector<double> covarBeta;
  std::vector<double> covarBetaOld;
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamBeta(const FixedData & fD) { init(fD); };

  virtual void setFill(std::vector<double> & probs,
		       const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD) const;

  virtual void updateFill(std::vector<double> & probs,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const;

};



#endif
