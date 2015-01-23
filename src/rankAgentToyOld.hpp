#ifndef RANK_AGENT_TOY_OLD_HPP__
#define RANK_AGENT_TOY_OLD_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include <limits>
#include "toyFeatures0.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "dataDepth.hpp"
#include "tuneParam.hpp"
#include "calcCentrality.hpp"


class RankToyOldTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  arma::colvec weights;

  int numChunks;
  int valReps;
};


template <class Model, class ModelParam>
class RankToyOldAgent : BaseAgent<Model,ModelParam> {
 public:
  RankToyOldAgent();
  
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const Model & m,
			ModelParam & mP);

  virtual void getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   const ModelParam & mP);

  std::vector<double> subgraph;
  
  arma::mat infFeat;
  arma::mat notFeat;
  arma::colvec infRanks;
  arma::colvec notRanks;
  
  static int numFeatures;

  int numAct;
  int numPre;

  RankToyOldTuneParam tp;

  std::string name;
};



#endif
