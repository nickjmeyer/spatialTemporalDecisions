#ifndef ANCHOR_MAN_HPP__
#define ANCHOR_MAN_HPP__



#include <vector>
#include <algorithm>
#include "m2NmOptim.hpp"


template <class System, class Agent>
class AnchorMan {
 public:
  AnchorMan();
  
  void addPar(const std::vector<double> & par);
  std::vector<std::vector<double> > parHist;

};



/* Make a switch optim that will cycle through the parHist as the "optimization"
   function so that existing runners can be used.  The runner to use would be
   the OptimRunner, but this class saves results to disk.  Re-write another
   OptimRunner that does not save any results.
*/




#endif
