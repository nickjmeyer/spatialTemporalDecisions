#ifndef AGENT_HPP__
#define AGENT_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"

template <class Model, class ModelParam>
class BaseAgent {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const Model & model,
			ModelParam & modelParam) = 0;

  std::string name;
};



int getNumPre(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD);

int getNumAct(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD);


#endif
