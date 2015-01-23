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
#include "toyFeatures0.hpp"


class RankTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  arma::colvec weights;

  int numChunks;
  int valReps;
};


template < class Features, class Model, class ModelParam>
class RankAgent : BaseAgent<Model,ModelParam> {
 public:
  RankAgent();
  
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

  arma::mat infFeat;
  arma::mat notFeat;
  arma::colvec infRanks;
  arma::colvec notRanks;
  
  static int numFeatures;

  int numAct;
  int numPre;

  RankTuneParam tp;

  std::string name;
};


#endif
