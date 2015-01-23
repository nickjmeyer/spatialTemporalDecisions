#ifndef M1_SIMPLE_OPTIM_HPP__
#define M1_SIMPLE_OPTIM_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"
#include "rankAgentToy.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "runner.hpp"

class M1SimpleOptimTunePar : public TuneParam{
 public:
  M1SimpleOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int mcReps; // number of Monte Carlo replications
};


template <class System, class Agent>
class M1SimpleOptim : BaseOptim<System,Agent>{
 public:
  M1SimpleOptim();
  
  virtual void optim(System system,
		     Agent & agent);
  M1SimpleOptimTunePar tp;

  std::string name;
};



#endif
