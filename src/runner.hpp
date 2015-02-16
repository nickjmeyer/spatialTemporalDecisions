#ifndef RUNNER_HPP__
#define RUNNER_HPP__


#include <omp.h>
#include <chrono>
#include "system.hpp"
#include "model.hpp"
#include "modelEbola.hpp"
#include "modelRange.hpp"
#include "agent.hpp"
#include "noTrtAgent.hpp"
#include "proximalAgent.hpp"
#include "myopicAgent.hpp"
#include "rankAgentToy.hpp"
#include "optim.hpp"
#include "m1SgdOptim.hpp"
#include "m2SaOptim.hpp"
#include "m1SimpleOptim.hpp"
#include "m1HybridOptim.hpp"



template <class S, class A>
class BaseRunner {
 public:
};




template <class S, class A>
class TrainRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     const int numReps, const int numPoints);
};




template <class S, class A>
class PlainRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     const int numReps, const int numPoints);

  std::pair<double,double> runEx(S system,
				 A agent,
				 const int numReps, const int numPoints);
};



template <class S, class A>
class VanillaRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     const int numReps, const int numPoints);
};


template <class S, class A>
class VanillaRunnerNS : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     const int numReps, const int numPoints);
};



template <class S, class A>
class FitOnlyRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     const int numReps, const int numPoints);
};



template <class S, class A, class Optim>
class OptimRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     Optim optim,
		     const int numReps, const int numPoints);
};



template <class S, class A, class Optim>
class OptimRunnerNS : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     Optim optim,
		     const int numReps, const int numPoints);
};



template <class S, class A, class Optim>
class TestRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     Optim optim,
		     const int numReps, const int numPoints);
};


template <class S, class A>
class TimerRunner : BaseRunner<S,A> {
 public:
  virtual double run(S system,
		     A agent,
		     const int numReps, const int numPoints);
};


#endif
