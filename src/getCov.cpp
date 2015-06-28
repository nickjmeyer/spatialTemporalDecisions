#include <iostream>
#include <cmath>
#include <cstdio>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>

extern "C" void getCov(double * const cov,
		       const double * const rnorm,
		       const double * const cenX,
		       const double * const cenY,
		       const int * const n_,
		       const double * rho_,
		       const double * tau_,
		       const double * eta_,
		       const int * p_,
		       const double * tol_);


void getSigmaSparse(Eigen::SparseMatrix<double> & sigma,
		    const double * const cenX,
		    const double * const cenY,
		    const int n,
		    const double rho,
		    const double tau,
		    const double eta,
		    const int p,
		    const double tol);


extern "C" void getCov(double * const cov,
		       const double * const rv,
		       const double * const cenX,
		       const double * const cenY,
		       const int * n_,
		       const double * rho_,
		       const double * tau_,
		       const double * eta_,
		       const int * p_,
		       const double * tol_){
  const int n = *n_;
  const double rho = *rho_;
  const double tau = *tau_;
  const double eta = *eta_;
  const int p = *p_;
  const double tol = *tol_;

  Eigen::SparseMatrix<double> sigma;
  getSigmaSparse(sigma,cenX,cenY,n,rho,tau,eta,p,tol);

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > llt;
  llt.compute(sigma);

  Eigen::Map<Eigen::VectorXd> covEig(cov,n*p);
  Eigen::Map<const Eigen::VectorXd> rvEig(rv,n*p);

  covEig = llt.matrixL() * rvEig;
}



void getSigmaSparse(Eigen::SparseMatrix<double> & sigma,
		    const double * const cenX,
		    const double * const cenY,
		    const int n,
		    const double rho,
		    const double tau,
		    const double eta,
		    const int p,
		    const double tol){
  int i,j,k,l,ind0,ind1;
  double val;
  double dist;

  sigma.resize(n*p,n*p);
  sigma.setZero();

  for(i = 0; i < n; ++i){
    for(j = 0; j < p; ++j){
      ind0 = i*p + j;
      for(k = 0; k < n; ++k){
	for(l = 0; l < p; ++l){
	  ind1 = k*p + l;
	  dist = std::sqrt(std::pow(cenX[i] - cenX[k],2.0) +
			   std::pow(cenY[i] - cenY[k],2.0));
	  val = rho*std::exp(-tau*dist - eta*std::abs(j-l));
	  if(val > tol){
	    sigma.insert(ind0,ind1) = val;
	  }
	}
      }
    }
  }
}
