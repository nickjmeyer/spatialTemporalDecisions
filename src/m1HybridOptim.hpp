#ifndef M1_HYBRID_OPTIM_HPP__
#define M1_HYBRID_OPTIM_HPP__


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

class M1HybridOptimTunePar : public TuneParam{
 public:
  M1HybridOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int topWeights;

  // for simple optim
  int mcRepsSimple;

  // for sgd optim
  int mcRepsSgd;
  int aSgd;
  int bSgd;
  double tolSgd;
};


template <class System, class Agent>
class M1HybridOptim : BaseOptim<System,Agent>{
 public:
  M1HybridOptim();
  
  virtual void optim(System system,
		     Agent & agent);
  M1HybridOptimTunePar tp;

  std::string name;
};



#endif
