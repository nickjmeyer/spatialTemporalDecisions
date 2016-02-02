#ifndef M2_Q_OPTIM_HPP__
#define M2_Q_OPTIM_HPP__


#include <vector>
#include <armadillo>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <gsl/gsl_bspline.h>
#include "timer.hpp"
#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "rankAgent.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"
#include "featuresInt.hpp"

#ifdef NJM_USE_MKL
#include "pardisoSymWrap.hpp"
#endif

class M2QEvalTunePar : public TuneParam{
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  int polReps; // reps to average over policy stochasticity
  int numNeigh; // number of Neighbors

  int dfLat; // bspline df for latitude
  int dfLong; // bspline df for longitude

  double gamma; // discount factor for bellman residual
  double lambda; // penalty coefficient

  int bootReps;
  double bootSize;

  double prune;

  double C,t,ell,muMin,A,B;  
};



template <class S, class A, class F, class M>
class M2QEval {
 public:

  M2QEval();

  double qFn(const SimData & sD, TrtData & tD,
	     const FixedData & fD, const DynamicData & dD,
	     M & m,
	     A a);
  
  double bellRes();

  // sets beta
  void solve();

  // set all pre-computed data
  void preCompData(const SimData & sD, const FixedData & fD);

  // fixed data for bellman residual
  // sets R,D0 (other stuff too, but this is the main stuff)
  void bellResFixData(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD,
		      M & m);

  // policy generated data for bellman residual
  // sets D1
  void bellResPolData(const int time,
		      const FixedData & fD,
		      M & m,
  		      A a);

  void buildRD();
  void buildRD(const std::vector<int> nodes);

  void buildD1();
  void buildD1(const std::vector<int> nodes);

  void getRD(Eigen::SparseVector<double> & R,
	     Eigen::SparseMatrix<double> & D);
  void setRD(const Eigen::SparseVector<double> & R,
	     const Eigen::SparseMatrix<double> & D);

  void getSolveDat(Eigen::VectorXd & mDtR,
		   Eigen::SparseMatrix<double> & DtD);
  void setSolveDat(const Eigen::VectorXd & mDtR,
		   const Eigen::SparseMatrix<double> & DtD);

  void fullRankPen();
  

  void tune();
  
  
  F f; // used to generate features
  
  std::vector<double> feat2Vec(const int numNodes,
			       const std::vector<int> & status);

  std::vector<Eigen::SparseMatrix<double>
	      >featToPsi(const std::vector<double> & feat);
  
  M2QEvalTunePar tp;

  int numFeat; // number of features w/ interactions  
  int lenPsi; // number of features w/ interactions AND neighbor average
  int dim; // (tp.dfLat * tp.dfLong + 1) * lenPsi;
  int numNodes;
  
  Eigen::VectorXd beta;

  // neighbors by distance
  std::vector<std::vector<int> > neighbors;

  // data for bellman residual
  std::vector<Eigen::SparseMatrix<double> > phiL;  // \lbrace \Phi_\ell \rbrace
  std::vector<std::vector<Eigen::SparseMatrix<double> > > phiPsiTL;

  std::vector<Eigen::SparseMatrix<double> > D0L;
  std::vector<Eigen::SparseMatrix<double> > D1L;
  std::vector<Eigen::SparseVector<double> > RL;
  
  Eigen::SparseVector<double> R;
  Eigen::VectorXd mDtR;
  Eigen::SparseMatrix<double> D0,D1,D,DtD,P,P2;
  std::vector<SimData> sD1T; // {sD}_{t=1}^T
  std::vector<DynamicData> dD1T; //{dD}_{t=1}^T
};



class M2QOptimTunePar : public TuneParam{
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  double C,t,ell,muMin,A,B;  
};




template <class S, class A, class F,
	  class M>
class M2QOptim : BaseOptim<S,A,M> {
 public:
  M2QOptim();
  
  void reset();

  virtual void optim(const S & system,
		     A & agent);
  
  M2QEval<System<M,M>,A,F,M> qEval;

  M2QOptimTunePar tp;
  
  std::string name;
};





#endif
