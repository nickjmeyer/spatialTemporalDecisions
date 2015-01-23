#ifndef M2_NM_OPTIM_OLD_HPP__
#define M2_NM_OPTIM_OLD_HPP__


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

class M2OldNmEvalTunePar : public TuneParam{
 public:
  M2OldNmEvalTunePar();
  
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




template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
class M2OldNmEval {
 public:
  double qFn(const SimData & sD, TrtData & tD,
	     const FixedData & fD, const DynamicData & dD,
	     const Model & m, ModelParam & mP,
	     Agent<Model,ModelParam> a);
  
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
  		      Agent<Model,ModelParam> a);
  
  // builds the spatial penalty
  Eigen::SparseMatrix<double> buildSpatialPen(const SimData & sD,
					      const FixedData & fD);
  Eigen::SparseMatrix<double> buildL2Pen(const int numNodes);

  std::vector<double> getFeatures(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  const Model & m,
				  ModelParam & mP);

  Eigen::SparseMatrix<double> featToPhi(const std::vector<double> & feat,
					const int numNodes);
  
  M2OldNmEvalTunePar tp;

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

  static const int numFeatures;
};



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
class M2OldNmOptim : BaseOptim<System,Agent,Model,ModelParam> {
 public:
  M2OldNmOptim();
  
  virtual void optim(System<Model,ModelParam> system,
		     Agent<Model,ModelParam> & agent);
  
  M2OldNmEval<System,Agent,Model,ModelParam> qEval;
  
  std::string name;
};



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
class M2OldNmData {
 public:
  M2OldNmEval<System,Agent,Model,ModelParam>  qEval;

  System<Model,ModelParam> s;
  Agent<Model,ModelParam> a;
};


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double M2OldNmObj(const gsl_vector * x, void * params){
  M2OldNmData<System,Agent,Model,ModelParam> * qD =
    static_cast<M2OldNmData<System,Agent,Model,ModelParam> *>(params);
  int i;
  for(i=0; i<qD->a.numFeatures; i++)
    qD->a.tp.weights(i) = gsl_vector_get(x,i);
  qD->qEval.bellResPolData(qD->s.sD.time,qD->s.fD,
			   qD->s.model,qD->s.estParam,
			   qD->a);
  qD->qEval.solve();
  return - qD->qEval.qFn(qD->s.sD,qD->s.tD,qD->s.fD,qD->s.dD,
			 qD->s.model,qD->s.estParam,qD->a);
}


#endif
