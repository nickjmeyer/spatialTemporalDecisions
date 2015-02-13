#ifndef OPTIM_HPP__
#define OPTIM_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgentToy.hpp"

template <class S, class A, class M, class MP>
class BaseOptim {
 public:
  virtual void optim(const S & system,
		     A & agent)=0;

  std::string name;
};



#endif
