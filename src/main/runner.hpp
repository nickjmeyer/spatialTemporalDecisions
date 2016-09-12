#ifndef RUNNER_HPP__
#define RUNNER_HPP__


#include <omp.h>
#include <chrono>
#include "runStats.hpp"
#include "starts.hpp"
#include "system.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "noTrtAgent.hpp"
#include "randomAgent.hpp"
#include "proximalAgent.hpp"
#include "myopicAgent.hpp"
#include "rankAgent.hpp"
#include "allAgent.hpp"
#include "optim.hpp"
#include "m1SpOptim.hpp"



template <class S, class A>
class BaseRunner {
 public:
};




template <class S, class A>
class TrainRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       const int numReps, const int numPoints);
};




template <class S, class A>
class PlainRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       const int numReps, const int numPoints);
};



template <class S, class A>
class VanillaRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       const int numReps, const int numPoints,
		       const Starts & starts);
};


template <class S, class A>
class VanillaRunnerNS : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       const int numReps, const int numPoints,
		       const Starts & starts);
};



template <class S, class A>
class FitOnlyRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       const int numReps, const int numPoints,
		       const Starts & starts);
};



template <class S, class A, class Optim>
class OptimRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       Optim optim,
		       const int numReps, const int numPoints,
		       const Starts & starts);
};



template <class S, class A, class Optim>
class OptimRunnerNS : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       Optim optim,
		       const int numReps, const int numPoints,
		       const Starts & starts);
};



template <class S, class A, class Optim>
class TuneRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       Optim optim,
		       const int numReps, const int numPoints);
};



template <class S, class A, class Optim>
class TestRunner : BaseRunner<S,A> {
 public:
  virtual RunStats run(S system,
		       A agent,
		       Optim optim,
		       const int numReps, const int numPoints);
};



#endif
