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
};

template <class System, class Agent>
class AnchorMan : BaseOptim<System,Agent> {
 public:
  AnchorMan();

  virtual void optim(System system,
		     Agent & agent);

  
};



/* Make a switch optim that will cycle through the parHist as the "optimization"
   function so that existing runners can be used.  The runner to use would be
   the OptimRunner, but this class saves results to disk.  Re-write another
   OptimRunner that does not save any results.
*/




#endif
