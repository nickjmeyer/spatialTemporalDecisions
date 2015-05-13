#ifndef PARAM_TRT_HPP__
#define PARAM_TRT_HPP__


#include "param.hpp"

class ParamTrt : public ParamBase {
 protected:
  std::vector<int> a;
  std::vector<int> p;
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamTrt( ) { };
  virtual ParamBase * clone() const {return new ParamTrt(*this);};

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
