#ifndef RANK_AGENT_HPP__
#define RANK_AGENT_HPP__


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
#include "calcCentrality.hpp"


class RankTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  arma::colvec weights;

  double jitterScale;
};


template < class F, class M, class MP>
class RankAgent : BaseAgent<M,MP> {
 public:
  RankAgent();

  void reset();
  
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

  RankTuneParam tp;

  std::string name;
};



#endif
