#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <vector>
// #include <armadillo>
// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
std::vector<double> getCovCpp(const std::vector<double> & rv,
			      const std::vector<double> & cenX,
			      const std::vector<double> & cenY,
			      const int n,
			      const double rho,
			      const double tau,
			      const double eta,
			      const int p);

// [[Rcpp::export]]
void getSigmaSparse(arma::mat & sigma,
		    const std::vector<double> & cenX,
		    const std::vector<double> & cenY,
		    const int n,
		    const double rho,
		    const double tau,
		    const double eta,
		    const int p);


// int main(){
//   std::ifstream ifs;
//   ifs.open("eta.txt");
//   double eta;
//   ifs >> eta;
//   ifs.close();

//   ifs.open("n.txt");
//   int n;
//   ifs >> n;
//   ifs.close();

//   n = 1000;

//   ifs.open("p.txt");
//   int p;
//   ifs >> p;
//   ifs.close();

//   ifs.open("rho.txt");
//   double rho;
//   ifs >> rho;
//   ifs.close();

//   ifs.open("rv.txt");
//   std::vector<double> rv;
//   for(int i = 0; i < n*p; ++i){
//     double val;
//     ifs >> val;
//     rv.push_back(val);
//   }
//   ifs.close();

//   ifs.open("tau.txt");
//   double tau;
//   ifs >> tau;
//   ifs.close();

//   ifs.open("x.txt");
//   std::vector<double> x;
//   for(int i = 0; i < n; ++i){
//     double val;
//     ifs >> val;
//     x.push_back(val);
//   }
//   ifs.close();

//   ifs.open("y.txt");
//   std::vector<double> y;
//   for(int i = 0; i < n; ++i){
//     double val;
//     ifs >> val;
//     y.push_back(val);
//   }
//   ifs.close();

//   std::vector<double> cov = getCovCpp(rv,x,y,n,rho,tau,eta,p);

//   return 0;
// }


std::vector<double> getCovCpp(const std::vector<double> & rv,
			      const std::vector<double> & cenX,
			      const std::vector<double> & cenY,
			      const int n,
			      const double rho,
			      const double tau,
			      const double eta,
			      const int p){
  arma::mat sigma;
  getSigmaSparse(sigma,cenX,cenY,n,rho,tau,eta,p);


  arma::colvec rvArma(rv);

  arma::mat r = arma::chol(sigma);

  return arma::conv_to<std::vector<double> >::from(r.t()*rvArma);
}



void getSigmaSparse(arma::mat & sigma,
		    const std::vector<double> & cenX,
		    const std::vector<double> & cenY,
		    const int n,
		    const double rho,
		    const double tau,
		    const double eta,
		    const int p){
  int i,j,k,l,ind0,ind1;
  double val;
  double dist;

  sigma.resize(n*p,n*p);
  sigma.zeros();

  for(i = 0; i < n; ++i){
    for(j = 0; j < p; ++j){
      ind0 = i*p + j;
      for(k = 0; k < n; ++k){
	for(l = 0; l < p; ++l){
	  ind1 = k*p + l;
	  dist = std::sqrt(std::pow(cenX.at(i) - cenX.at(k),2.0) +
			   std::pow(cenY.at(i) - cenY.at(k),2.0));
	  val = rho*std::exp(-tau*dist - eta*std::abs(j-l));

	  sigma(ind0,ind1) = val;
	}
      }
    }
  }
}
