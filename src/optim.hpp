#ifndef OPTIM_HPP__
#define OPTIM_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"

template <class System, class Agent>
class BaseOptim {
 public:
  virtual void optim(System system,
		     Agent & agent)=0;

  std::string name;
};



#endif
