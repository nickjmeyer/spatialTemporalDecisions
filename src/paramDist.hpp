#ifndef PARAM_DIST_HPP__
#define PARAM_DIST_HPP__

#include "param.hpp"

class ParamDist : public ParamBase {
 protected:
  std::vector<double> dist;
  std::vector<double> alphaDist;
  int numNodes;

  virtual unsigned int initParsSize(const FixedData & fD);

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamDist() { };
  virtual ParamBase * clone() const {return new ParamDist(*this);};

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
