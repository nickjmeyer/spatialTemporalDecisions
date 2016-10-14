#ifndef OPTIM_HPP
#define OPTIM_HPP


#include "data.hpp"
#include "model.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"

template <class S, class A, class M>
class BaseOptim {
 public:
  virtual void optim(const S & system,
		     A & agent)=0;

  std::string name;
};



#endif
