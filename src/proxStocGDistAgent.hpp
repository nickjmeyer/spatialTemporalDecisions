#ifndef PROXIMAL_STOCHASTIC_G_DIST_AGENT_HPP__
#define PROXIMAL_STOCHASTIC_G_DIST_AGENT_HPP__


#include <vector>
#include <queue>
#include <limits>
#include "data.hpp"
#include "model.hpp"
#include "modelGravityGDist.hpp"
#include "agent.hpp"
#include "tuneParam.hpp"
#include "features.hpp"
#include "runStats.hpp"

class ProxStocGDistTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  double corrGoal;
};
  

template <class M>
class ProxStocGDistAgent : public BaseAgent<M> {
 public:
  ProxStocGDistAgent();
  
  virtual void applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model);

  virtual void reset();

  int numAct;
  int numPre;

  static std::string name;

  ProxStocGDistTuneParam tp;
};


double optCorr(const std::vector<double> vals,
	       const double goal);

std::vector<double> getScores(std::vector<double> vals);

std::vector<std::vector<double> >
appScale(std::vector<std::vector<double> > vals,
	 const double scale);

double getCorr(const std::vector<double> vals,
	       const std::vector<std::vector<double> > jitter);

double getCorr(const std::vector<double> vals0,
	       const std::vector<double> vals1);

std::vector<double> getRank(const std::vector<double> scores);



#endif
