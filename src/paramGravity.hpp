#ifndef PARAM_GRAVITY_HPP__
#define PARAM_GRAVITY_HPP__


#include "param.hpp"

class ParamGravity : public ParamBase {
 protected:
  std::vector<double> parsOld;
  std::vector<double> grav;
  std::vector<double> gravOld;
  std::vector<double> dist;
  std::vector<double> cc;
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamGravity(const FixedData & fD) { init(fD); };

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
