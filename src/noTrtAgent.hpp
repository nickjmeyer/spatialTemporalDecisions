#ifndef NO_TRT_AGENT_HPP__
#define NO_TRT_AGENT_HPP__

#include "agent.hpp"
#include "features.hpp"


template <class M, class MP>
class NoTrt : public BaseAgent<M,MP> {
 public:
  void applyTrt(const SimData & sD,
		TrtData & tD,
		const FixedData & fD,
		const DynamicData & dD,
		const M & model,
		MP & modelParam);

  static const std::string name;
};


#endif
