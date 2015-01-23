#ifndef Q_OPTIM_GSL_HPP__
#define Q_OPTIM_GSL_HPP__


#include <vector>
#include <armadillo>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "model.hpp"
#include "modelParam.hpp"
#include "system.hpp"
#include "agent.hpp"
#include "rankAgent.hpp"
#include "rankAgentToy.hpp"
#include "optim.hpp"
#include "tuneParam.hpp"

class QOptimTunePar : public TuneParam{
 public:
  QOptimTunePar();
  
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  double radius; // radius for neighbors
  int polReps; // reps to average over policy stochasticity
  
  double jitter; // standard deviation of Gaussian noise
  double tol; // convergence tolerance
  double rate; // learning rate
  double rateDecay; // learning rate decay
};


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
class QOptim : BaseOptim<System,Agent,Model,ModelParam> {
 public:
  virtual void optim(System<Model,ModelParam> system,
		     Agent<Model,ModelParam> & agent);

  double Qfn(System<Model,ModelParam> s,
	     Agent<Model,ModelParam> a);

  double bmRes();
  arma::colvec bmResG(); // G means gradient


  double spPen();
  arma::colvec spPenG(); // G means gradient


  // used for bellman residual
  void setRD(const SimData & sD,
	     const TrtData & tD,
	     const FixedData & fD,
	     const DynamicData & dD,
	     const Model & m,
	     ModelParam & mP,
	     Agent<Model,ModelParam> & agent);
  
  // used for spatial penalty
  void setHMG(const SimData & sD, const FixedData & fD);

  arma::colvec getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);
  
  QOptimTunePar tp;

  arma::colvec beta;

  // for bellman residual
  arma::colvec R;
  arma::mat D;
  
  // \sum_{\ell} (H_\ell - G_\ell)^\T(H_\ell - G_\ell)
  arma::sp_mat hmg;

  static const int numFeatures;
};


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double QOptimObjFn(const gsl_vector * x, void * param);

template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptimObjG(const gsl_vector * x, void * param, gsl_vector * g);

template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptimObjFnG(const gsl_vector * x, void * param,
		  double * f, gsl_vector * g);


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
class QOptimData {
 public:
  QOptim<System,Agent,Model,ModelParam> q;
  System<Model,ModelParam> s;
  Agent<Model,ModelParam> a;
};


#endif
