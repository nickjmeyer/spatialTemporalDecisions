#ifndef PARAM_RADIUS_HPP__
#define PARAM_RADIUS_HPP__

#include "param.hpp"

class ParamRadius : public ParamBase {
 protected:
  std::vector<double> logDist;
  std::vector<int> beyond;
  int numNodes;

  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual std::vector<bool> initToScale();

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamRadius() { };
  virtual ParamBase * clone() const {return new ParamRadius(*this);};

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
