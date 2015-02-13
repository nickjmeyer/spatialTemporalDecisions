#ifndef MYOPIC_AGENT_HPP__
#define MYOPIC_AGENT_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "features.hpp"
#include "modelParamRange.hpp"
#include "modelParamCave.hpp"
#include "modelRange.hpp"
#include "modelCave.hpp"

template <class M, class MP>
class MyopicAgent : BaseAgent<M,MP> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const M & model,
			MP & modelParam);

  arma::colvec infFeat;
  arma::colvec notFeat;

  int numAct;
  int numPre;

  static std::string name;
};


#endif
