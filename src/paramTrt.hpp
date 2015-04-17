#ifndef PARAM_TRT_HPP__
#define PARAM_TRT_HPP__


#include "param.hpp"

class ParamTrt : public ParamBase {
 protected:
  std::vector<double> parsOld;
  std::vector<int> a;
  std::vector<int> p;
  std::vector<std::pair<int,double> > trt;
  std::vector<std::pair<int,double> > trtOld;
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamTrt(const FixedData & fD) { init(fD); };

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
