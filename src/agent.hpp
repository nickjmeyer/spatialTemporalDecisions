#ifndef AGENT_HPP__
#define AGENT_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"

template <class M, class MP>
class BaseAgent {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const M & model,
			MP & modelParam) = 0;

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
