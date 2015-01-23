#ifndef M1_NM_OPTIM_HPP__
#define M1_NM_OPTIM_HPP__


#include <gsl/gsl_multimin.h>
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

class M1NmOptimTunePar : public TuneParam{
 public:
  M1NmOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int mcReps; // number of Monte Carlo replications
  double tol; // convergence tolerance
  double ss; // starting step size
};


template <class System,class Agent>
class M1NmOptim : BaseOptim<System,Agent>{
 public:
  M1NmOptim();
  virtual void optim(System system,
		     Agent & agent);

  M1NmOptimTunePar tp;
  
  std::string name;
};




template <class System, class Agent>
class M1NmData {
 public:
  System s;
  
  Agent a;
  
  PlainRunner<System,Agent> r;

  int numReps;
  int numYears;
};



template <class System, class Agent>
double M1NmObj(const gsl_vector * x, void * params){
  M1NmData<System,Agent> * d =
    static_cast<M1NmData<System,Agent> *>(params);
  int i;
  for(i=0; i<d->a.f.numFeatures; i++)
    d->a.tp.weights(i) = gsl_vector_get(x,i);
  fixRandomSeed(1);
  double val = d->r.run(d->s,d->a,d->numReps,d->numYears);
  fixRandomSeed(0);
  return val;
}




#endif
