#ifndef PARAM_INTERCEPT_HPP__
#define PARAM_INTERCEPT_HPP__

#include "param.hpp"

class ParamIntercept : public ParamBase {
 protected:
  double parsOld;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamIntercept(const FixedData & fD) { init(fD); };

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
