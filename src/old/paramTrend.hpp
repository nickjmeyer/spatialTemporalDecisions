#ifndef PARAM_TREND_HPP__
#define PARAM_TREND_HPP__

#include "param.hpp"

class ParamTrend : public ParamBase {
 protected:
  unsigned int prevTime;
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);

  virtual std::vector<std::string> initNames();
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamTrend() { };

  virtual ParamBase * clone() const {return new ParamTrend(*this);};

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
