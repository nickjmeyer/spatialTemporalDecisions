#ifndef MYOPIC_AGENT_HPP__
#define MYOPIC_AGENT_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "features.hpp"

template <class Model, class ModelParam>
class MyopicAgent : BaseAgent<Model,ModelParam> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const Model & model,
			ModelParam & modelParam);

  arma::colvec infFeat;
  arma::colvec notFeat;

  int numAct;
  int numPre;

  static const std::string name;
};


#endif
