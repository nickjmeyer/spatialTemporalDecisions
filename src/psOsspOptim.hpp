#ifndef PS_OSS_OPTIM_HPP__
#define PS_OSS_OPTIM_HPP__


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

class PsOsspOptimTunePar : public TuneParam{
 public:
  PsOsspOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int N,B,mcReps;
  double corrGoal;
};


template <class S, class A, class M>
class PsOsspOptim : BaseOptim<S,A,M>{
 public:
  PsOsspOptim();

  void reset();
  
  virtual void optim(const S & system,
		     A & agent);
  
  PsOsspOptimTunePar tp;

  std::string name;
};



#endif
