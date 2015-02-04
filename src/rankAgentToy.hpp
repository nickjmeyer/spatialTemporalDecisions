#ifndef RANK_AGENT_TOY_HPP__
#define RANK_AGENT_TOY_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include <limits>
#include "system.hpp"
#include "agent.hpp"
#include "dataDepth.hpp"
#include "tuneParam.hpp"
#include "features.hpp"
#include "toyFeatures0.hpp"
#include "toyFeatures1.hpp"
#include "toyFeatures2.hpp"
#include "calcCentrality.hpp"


class RankToyTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  arma::colvec weights;

  int numChunks;
};


template < class F, class M, class MP>
class RankToyAgent : BaseAgent<M,MP> {
 public:
  RankToyAgent();
  
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const M & m,
			MP & mP);

  F f;

  arma::colvec infRanks;
  arma::colvec notRanks;
  
  int numAct;
  int numPre;

  RankToyTuneParam tp;

  std::string name;
};



#endif
