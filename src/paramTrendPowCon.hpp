#ifndef PARAM_TREND_POW_CON_HPP__
#define PARAM_TREND_POW_CON_HPP__

#include "param.hpp"

class ParamTrendPowCon : public ParamBase {
 protected:
  unsigned int prevTime;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();  

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamTrendPowCon() { };

  virtual ParamBase * clone() const {return new ParamTrendPowCon(*this);};

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
  
  virtual std::vector<double> partial2(const int notNode,
				       const int infNode,
				       const SimData & sD,
				       const TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD);
  
};





#endif
