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
#include "rankAgent.hpp"
#include "rankAgentToy.hpp"
#include "optim.hpp"
#include "m1SgdOptim.hpp"
#include "m2NmOptim.hpp"
#include "m1SimpleOptim.hpp"
#include "m1HybridOptim.hpp"



template <class System, class Agent>
class BaseRunner {
 public:
};




template <class System, class Agent>
class TrainRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     const int numReps, const int numPoints);
};




template <class System, class Agent>
class PlainRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     const int numReps, const int numPoints);

  std::pair<double,double> runEx(System system,
				 Agent agent,
				 const int numReps, const int numPoints);
};



template <class System, class Agent>
class VanillaRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     const int numReps, const int numPoints);
};



template <class System, class Agent>
class FitOnlyRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     const int numReps, const int numPoints);
};



template <class System, class Agent, class Optim>
class OptimRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     Optim optim,
		     const int numReps, const int numPoints);
};



template <class System, class Agent, class Optim>
class OptimRunnerNS : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     Optim optim,
		     const int numReps, const int numPoints);
};



template <class System, class Agent, class Optim>
class TestRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     Optim optim,
		     const int numReps, const int numPoints);
};


template <class System, class Agent>
class TimerRunner : BaseRunner<System,Agent> {
 public:
  virtual double run(System system,
		     Agent agent,
		     const int numReps, const int numPoints);
};


#endif
