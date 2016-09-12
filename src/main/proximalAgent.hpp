#ifndef PROXIMAL_AGENT_HPP__
#define PROXIMAL_AGENT_HPP__


#include <vector>
#include <queue>
#include <limits>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "features.hpp"


template <class M>
class ProximalAgent : public BaseAgent<M> {
 public:
  ProximalAgent();

  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model);

  int numAct;
  int numPre;

  virtual void setEdgeToEdge(const bool edgeToEdge);

  bool edgeToEdge;

  static std::string name;
};



#endif
