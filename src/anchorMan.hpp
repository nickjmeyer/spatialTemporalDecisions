#ifndef ANCHOR_MAN_HPP__
#define ANCHOR_MAN_HPP__



#include <vector>
#include <algorithm>
#include "m2NmOptim.hpp"


class AnchorManTunePar : public TuneParam {
 public:
  AnchorManTunePar();

  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int numSamples;
  double cutoff;
  int freq;
};

template <class System, class Agent, class Features,
	  class Model, class ModelParam>
class AnchorMan : BaseOptim<System,Agent> {
 public:
  AnchorMan();

  M1HybridOptim<Sytem,Agent> m1Opt;
  M2NmOptim<system,Agent,Features,Model,ModelParam> m2Opt;

  std::vector<double> m1W,m2W;

  virtual void optim(System system,
		     Agent & agent);

  int toSwitch(System system,
	       Agent & agent);
  
  double sampleNull(System system,
		    Agent & agent,
		    const int numYears);

  int switched;
};



#endif
