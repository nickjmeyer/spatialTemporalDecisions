#ifndef NO_TRT_AGENT_HPP
#define NO_TRT_AGENT_HPP

#include "agent.hpp"
#include "features.hpp"


template <class M>
class NoTrt : public BaseAgent<M> {
 public:
  void applyTrt(const SimData & sD,
		TrtData & tD,
		const FixedData & fD,
		const DynamicData & dD,
		M & model);

  static std::string name;
};


#endif
