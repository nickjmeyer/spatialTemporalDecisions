#ifndef OSSP_AGENT_HPP__
#define OSSP_AGENT_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include <limits>
#include "system.hpp"
#include "agent.hpp"
#include "dataDepth.hpp"
#include "tuneParam.hpp"
#include "features.hpp"
#include "toyFeatures2.hpp"


class OsspAgentTuneParam : public TuneParam {
 public:
  OsspAgentTuneParam();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  double lambda;
};


template <class M>
class OsspAgent : public BaseAgent<M> {
 public:
  OsspAgent();

  void reset();

  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & m);

  std::vector<std::vector<int> > pCand;
  std::vector<std::vector<int> > aCand;

  std::vector<double> qvalues;

  OsspAgentTuneParam tp;

  std::string name;
};




#endif
