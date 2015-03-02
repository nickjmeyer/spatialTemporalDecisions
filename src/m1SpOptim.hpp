#ifndef M1_SP_OPTIM_HPP__
#define M1_SP_OPTIM_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "runner.hpp"

class M1SpOptimTunePar : public TuneParam{
 public:
  M1SpOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int mcReps;

  double C,t,ell,muMin,A,B;

  int tune;
};


template <class S, class A, class M, class MP>
class M1SpOptim : BaseOptim<S,A,M,MP>{
 public:
  M1SpOptim();

  void reset();
  
  virtual void optim(const S & system,
		     A & agent);
  virtual void tune(const System<M,MP,M,MP> & system,
		    A agent);
  
  M1SpOptimTunePar tp;

  std::string name;
};



#endif
