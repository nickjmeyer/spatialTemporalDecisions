#ifndef PARAM_GRAVITY_G_DIST_HPP__
#define PARAM_GRAVITY_G_DIST_HPP__


#include "param.hpp"

class ParamGravityGDist : public ParamBase {
 protected:
  std::vector<double> grav;
  std::vector<double> dist;
  std::vector<double> cc;
  int numNodes;
  
  virtual unsigned int initParsSize(const FixedData & fD);

  virtual std::vector<std::string> initNames();

  virtual std::vector<bool> initToScale();

  virtual void initInternal(const FixedData & fD);
  
  virtual void updateBefore();
  
  virtual void updateAfter();

 public:
  ParamGravityGDist() { };
  virtual ParamBase * clone() const {return new ParamGravityGDist(*this);};

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