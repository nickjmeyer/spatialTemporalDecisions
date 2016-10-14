#ifndef ALL_AGENT_HPP
#define ALL_AGENT_HPP


#include <armadillo>
#include <vector>
#include <queue>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "features.hpp"

template <class M>
class AllAgent : public BaseAgent<M> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model);

  static std::string name;
};


#endif
