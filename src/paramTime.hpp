#ifndef PARAM_TIME_HPP__
#define PARAM_TIME_HPP__


#include "param.hpp"

class ParamTime : public ParamBase {
 protected:
  std::vector<int> time;
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamTime( ) { };
  virtual ParamBase * clone() const {return new ParamTime(*this);};

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
