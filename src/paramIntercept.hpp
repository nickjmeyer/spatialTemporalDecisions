#ifndef PARAM_INTERCEPT_HPP__
#define PARAM_INTERCEPT_HPP__

#include "param.hpp"

class ParamIntercept : public ParamBase {
 protected:
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamIntercept() { };

  virtual ParamBase * clone() const {return new ParamIntercept(*this);};

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
};





#endif
