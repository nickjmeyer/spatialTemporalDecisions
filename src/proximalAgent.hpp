#ifndef PROXIMAL_AGENT_HPP__
#define PROXIMAL_AGENT_HPP__


#include <vector>
#include <queue>
#include <limits>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "features.hpp"
#include "modelParamRange.hpp"
#include "modelRange.hpp"


template <class M, class MP>
class ProximalAgent : BaseAgent<M, MP> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const M & model,
			MP & modelParam);

  int numAct;
  int numPre;

  static std::string name;
};



#endif
