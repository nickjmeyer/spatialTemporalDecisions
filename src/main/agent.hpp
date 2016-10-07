#ifndef AGENT_HPP__
#define AGENT_HPP__


#include "data.hpp"
#include "model.hpp"

template <class M>
class BaseAgent {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model) = 0;

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


void checkForValidTrt(
  const SimData & sD, const TrtData & tD,
    const FixedData & fD, const DynamicData & dD);



#endif
