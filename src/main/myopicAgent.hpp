#ifndef MYOPIC_AGENT_HPP__
#define MYOPIC_AGENT_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include "data.hpp"
#include "model.hpp"
#include "agent.hpp"
#include "features.hpp"

template <class M>
class MyopicAgent : public BaseAgent<M> {
 public:
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model);

  arma::colvec infFeat;
  arma::colvec notFeat;

  int numAct;
  int numPre;

  static std::string name;
};


#endif
