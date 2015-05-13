#ifndef PARAM_TIME_EXP_CAVES_HPP__
#define PARAM_TIME_EXP_CAVES_HPP__


#include "param.hpp"

class ParamTimeExpCaves : public ParamBase {
 protected:
  std::vector<int> time;
  std::vector<double> iPropCaves; // (max(caves)+1)/(caves + 1)
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamTimeExpCaves( ) { };
  virtual ParamBase * clone() const {return new ParamTimeExpCaves(*this);};

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
