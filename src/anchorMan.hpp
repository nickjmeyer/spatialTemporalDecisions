#ifndef ANCHOR_MAN_HPP__
#define ANCHOR_MAN_HPP__



#include <vector>
#include <algorithm>
#include "utilities.hpp"
#include "data.hpp"
#include "system.hpp"
#include "model.hpp"
#include "modelRange.hpp"
#include "modelParam.hpp"
#include "modelParamRange.hpp"
#include "m1HybridOptim.hpp"
#include "m2NmOptim.hpp"
#include "rankAgentToy.hpp"


class AnchorManTunePar : public TuneParam {
 public:
  AnchorManTunePar();

  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int numSamples;
  double cutoff;
  int freq;
};


template <class S, class A, class F,
	  class M, class MP>
class AnchorMan : BaseOptim<S,A,M,MP> {
 public:
  AnchorMan();

  M1HybridOptim<System<M,MP,M,MP>,A,M,MP> m1Opt;
  M2NmOptim<System<M,MP,M,MP>,A,F,M,MP> m2Opt;

  std::vector<double> m1W,m2W;

  virtual void optim(const S & system,
		     A & agent);

  int toSwitch(System<M,MP,M,MP> & system,
	       A & agent, const int T);
  
  AnchorManTunePar tp;

  int switched;

  std::string name;
};



#endif
