#ifndef M1_SGD_OPTIM_HPP__
#define M1_SGD_OPTIM_HPP__


#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgentToy.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "runner.hpp"

class M1SgdOptimTunePar : public TuneParam{
 public:
  M1SgdOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  double jitter; // standard deviation of Gaussian noise
  int mcReps; // number of Monte Carlo replications
  double tol; // convergence tolerance
  double rate; // learning rate
  double rateDecay; // learning rate decay

  double momRate,a,b;
  
  int tune;

};


template <class S, class A, class M, class MP>
class M1SgdOptim : BaseOptim<S,A,M,MP>{
 public:
  M1SgdOptim();
  virtual void optim(const S & system,
		     A & agent);
  virtual void tune(const System<M,MP,M,MP> & system,
		    A agent);
  
  M1SgdOptimTunePar tp;

  std::string name;
};



#endif
