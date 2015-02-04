#ifndef RANDOM_AGENT_HPP__
#define RANDOM_AGENT_HPP__


#include <vector>
#include <queue>
#include <limits>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"



template <class M, class MP>
class RandomAgent : BaseAgent<M, MP> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const M & model,
			MP & modelParam);

  int numAct;
  int numPre;

  static const std::string name;
};



#endif
