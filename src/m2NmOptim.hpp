#ifndef M2_NM_OPTIM_HPP__
#define M2_NM_OPTIM_HPP__


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

class M2NmEvalTunePar : public TuneParam{
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
};



template <class System, class Agent, class Features,
	  class Model, class ModelParam>
class M2NmEval {
 public:

  M2NmEval();

  double qFn(const SimData & sD, TrtData & tD,
	     const FixedData & fD, const DynamicData & dD,
	     const Model & m, ModelParam & mP,
	     Agent a);
  
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
		      const Model & m,
		      ModelParam & mP);

  // policy generated data for bellman residual
  // sets D1
  void bellResPolData(const int time,
		      const FixedData & fD,
		      const Model & m,
  		      ModelParam & mP,
  		      Agent a);
  
  // builds the spatial penalty
  Eigen::SparseMatrix<double> buildSpatialPen(const SimData & sD,
					      const FixedData & fD);
  Eigen::SparseMatrix<double> buildL2Pen(const int numNodes);

  Features f; // used to generate features
  std::vector<double> feat2Vec(const int numNodes,
			       const std::vector<int> & status);

  Eigen::SparseMatrix<double> featToPhi(const std::vector<double> & feat,
					const int numNodes);
  
  M2NmEvalTunePar tp;

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



template <class System, class Agent, class Features,
	  class Model,class ModelParam>
class M2NmOptim : BaseOptim<System,Agent> {
 public:
  M2NmOptim();
  
  virtual void optim(System system,
		     Agent & agent);
  
  M2NmEval<System,Agent,Features,Model,ModelParam> qEval;
  
  std::string name;
};



template <class System, class Agent, class Features,
	  class Model,class ModelParam>
class M2NmData {
 public:
  M2NmEval<System,Agent,Features,Model,ModelParam>  qEval;

  System s;
  Agent a;
};


template <class System, class Agent, class Features,
	  class Model,class ModelParam>
double M2NmObj(const gsl_vector * x, void * params){
  M2NmData<System,Agent,Features,Model,ModelParam> * qD =
    static_cast<M2NmData<System,Agent,Features,Model,ModelParam> *>(params);
  int i;
  for(i=0; i<qD->a.f.numFeatures; i++)
    qD->a.tp.weights(i) = gsl_vector_get(x,i);
  qD->qEval.bellResPolData(qD->s.sD.time,qD->s.fD,
			   qD->s.model,qD->s.estParam,
			   qD->a);
  std::cout << std::endl << qD->qEval.D1.sum() << std::endl;
  qD->qEval.solve();
  return - qD->qEval.qFn(qD->s.sD,qD->s.tD,qD->s.fD,qD->s.dD,
			 qD->s.model,qD->s.estParam,qD->a);
}


#endif
