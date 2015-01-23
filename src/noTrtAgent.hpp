#ifndef NO_TRT_AGENT_HPP__
#define NO_TRT_AGENT_HPP__

#include "agent.hpp"
#include "features.hpp"


template <class Model, class ModelParam>
class NoTrt : public BaseAgent<Model,ModelParam> {
 public:
  void applyTrt(const SimData & sD,
		TrtData & tD,
		const FixedData & fD,
		const DynamicData & dD,
		const Model & model,
		ModelParam & modelParam);

  static const std::string name;
};


#endif
