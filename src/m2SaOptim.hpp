#ifndef M2_SA_OPTIM_HPP__
#define M2_SA_OPTIM_HPP__


#include <vector>
#include <armadillo>
#include <gsl/gsl_multimin.h>
#include <eigen3/Eigen/Sparse>
#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "proximalAgent.hpp"
#include "rankAgent.hpp"
#include "rankAgentToy.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "featuresInt.hpp"

class M2SaEvalTunePar : public TuneParam{
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int polReps; // reps to average over policy stochasticity
  int valReps;
  int numNeigh; // number of Neighbors

  double gamma; // discount factor for bellman residual
  double lambda; // penalty coefficient
  
  double jitter; // standard deviation of Gaussian noise
  double tol; // convergence tolerance
  double rate; // learning rate
  double rateDecay; // learning rate decay

  double sdStart;
  double sdStop;
  double sdDecay;
  double sdJump;
};



template <class S, class A, class F,
	  class M, class MP>
class M2SaEval {
 public:

  M2SaEval();

  double qFn(const SimData & sD, TrtData & tD,
	     const FixedData & fD, const DynamicData & dD,
	     const M & m, MP & mP,
	     A a);
  
  double bellRes();

  // sets beta
  void solve();

  // set all pre-computed data
  void preCompData(const SimData & sD, const FixedData & fD);

  // fixed data for bellman residual
  // sets R,D0,deltaQ (other stuff too, but this is the main stuff)
  void bellResFixData(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      const M & m,
		      MP & mP);

  // policy generated data for bellman residual
  // sets D1
  void bellResPolData(const int time,
		      const FixedData & fD,
		      const M & m,
  		      MP & mP,
  		      A a);
  
  // builds the spatial penalty
  Eigen::SparseMatrix<double> buildSpatialPen(const SimData & sD,
					      const FixedData & fD);
  Eigen::SparseMatrix<double> buildL2Pen(const int numNodes);

  F f; // used to generate features
  std::vector<double> feat2Vec(const int numNodes,
			       const std::vector<int> & status);

  Eigen::SparseMatrix<double> featToPhi(const std::vector<double> & feat,
					const int numNodes);
  
  M2SaEvalTunePar tp;

  int K; // number of features w/ interactions but w/ neighbor average
  
  Eigen::VectorXd beta;

  // spatial pentaly + L2 penalty
  Eigen::SparseMatrix<double> P;

  // neighbors by distance
  std::vector<std::vector<int> > neighbors;

  // data for bellman residual
  std::vector<Eigen::SparseMatrix<double> > deltaQ;
  Eigen::VectorXd R;
  Eigen::SparseMatrix<double> D0,D1,D;
  std::vector<SimData> sD1T; // {sD}_{t=1}^T
  std::vector<DynamicData> dD1T; //{dD}_{t=1}^T
};



template <class S, class A, class F,
	  class M,class MP>
class M2SaOptim : BaseOptim<S,A,M,MP> {
 public:
  M2SaOptim();
  
  virtual void optim(const S & system,
		     A & agent);
  
  M2SaEval<System<M,MP,M,MP>,A,F,M,MP> qEval;
  
  std::string name;
};





#endif
