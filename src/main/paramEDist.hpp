#ifndef PARAM_E_DIST_HPP
#define PARAM_E_DIST_HPP

#include "param.hpp"

class ParamEDist : public ParamBase {
 protected:
  std::vector<double> dist;
  std::vector<double> alphaDist;
  int numNodes;

  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual void initInternal(const FixedData & fD);

  virtual void updateBefore();

  virtual void updateAfter();

 public:
  ParamEDist() { };
  virtual ParamBase * clone() const {return new ParamEDist(*this);};

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
