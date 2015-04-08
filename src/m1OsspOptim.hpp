#ifndef M1_OSS_OPTIM_HPP__
#define M1_OSS_OPTIM_HPP__


#include <set>
#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "runner.hpp"

class M1OsspOptimTunePar : public TuneParam{
 public:
  M1OsspOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int N;
};


template <class S, class A, class M>
class M1OsspOptim : BaseOptim<S,A,M>{
 public:
  M1OsspOptim();

  void reset();
  
  virtual void optim(const S & system,
		     A & agent);
  
  virtual void tune(const System<M,M> & system,
		    A agent);
  
  M1OsspOptimTunePar tp;

  std::string name;
};



#endif
