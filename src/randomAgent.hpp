#ifndef RANDOM_AGENT_HPP__
#define RANDOM_AGENT_HPP__


#include <vector>
#include <queue>
#include <limits>
#include "data.hpp"
#include "model.hpp"
#include "modelGravity.hpp"
#include "agent.hpp"



template <class M>
class RandomAgent : BaseAgent<M> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model);

  int numAct;
  int numPre;

  static const std::string name;
};



#endif
