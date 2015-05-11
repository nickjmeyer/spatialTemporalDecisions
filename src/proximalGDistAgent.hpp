#ifndef PROXIMAL_G_DIST_AGENT_HPP__
#define PROXIMAL_G_DIST_AGENT_HPP__


#include <vector>
#include <queue>
#include <limits>
#include "data.hpp"
#include "model.hpp"
#include "modelGravityGDist.hpp"
#include "agent.hpp"
#include "features.hpp"


template <class M>
class ProximalGDistAgent : public BaseAgent<M> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model);

  int numAct;
  int numPre;

  static std::string name;
};



#endif
