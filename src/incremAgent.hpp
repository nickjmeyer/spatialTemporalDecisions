#ifndef INCREM_AGENT_HPP__
#define INCREM_AGENT_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include <limits>
#include "system.hpp"
#include "timer.hpp"
#include "agent.hpp"
#include "tuneParam.hpp"
#include "runner.hpp"
#include "proximalGDistAgent.hpp"
#include "nullOptim.hpp"


class IncremTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int N;
  int mcReps;
};


template <class M, class A, class O>
class IncremAgent : public BaseAgent<M> {
 public:
  IncremAgent();

  void reset();
  
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & m);

  virtual double eval(System<M,M> s,
		      A a);

  int numAct;
  int numPre;

  IncremTuneParam tp;

  std::string name;
};



#endif
