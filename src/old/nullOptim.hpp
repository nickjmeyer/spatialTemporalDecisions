#ifndef NULL_OPTIM_HPP__
#define NULL_OPTIM_HPP__


#include <set>
#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"
#include "osspAgent.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "runner.hpp"
#include "proxStocGDistAgent.hpp"


template <class S, class A, class M>
class NullOptim : BaseOptim<S,A,M>{
 public:
  NullOptim();

  void reset();
  
  virtual void optim(const S & system,
		     A & agent);
  
  std::string name;
};



#endif
