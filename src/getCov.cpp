#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <vector>
// #include <Rcpp.h>

// [[Rcpp::export]]
std::vector<double> getCovCpp(const std::vector<double> & rv,
			      const std::vector<double> & cenX,
			      const std::vector<double> & cenY,
			      const int n,
			      const double rho,
			      const double tau,
			      const double eta,
			      const int p,
			      const double tol);


void getSigmaSparse(Eigen::SparseMatrix<double> & sigma,
		    const std::vector<double> & cenX,
		    const std::vector<double> & cenY,
		    const int n,
		    const double rho,
		    const double tau,
		    const double eta,
		    const int p,
		    const double tol);


int main(){
  std::ifstream ifs;
  ifs.open("eta.txt");
  double eta;
  ifs >> eta;
  ifs.close();

  ifs.open("n.txt");
  int n;
  ifs >> n;
  ifs.close();

  ifs.open("p.txt");
  int p;
  ifs >> p;
  ifs.close();

  ifs.open("rho.txt");
  double rho;
  ifs >> rho;
  ifs.close();

  ifs.open("rv.txt");
  std::vector<double> rv;
  while(ifs.good()){
    double val;
    ifs >> val;
    rv.push_back(val);
  }
  ifs.close();

  ifs.open("tau.txt");
  double tau;
  ifs >> tau;
  ifs.close();

  ifs.open("tol.txt");
  double tol;
  ifs >> tol;
  ifs.close();

  ifs.open("x.txt");
  std::vector<double> x;
  while(ifs.good()){
    double val;
    ifs >> val;
    x.push_back(val);
  }
  ifs.close();

  ifs.open("y.txt");
  std::vector<double> y;
  while(ifs.good()){
    double val;
    ifs >> val;
    y.push_back(val);
  }
  ifs.close();

  std::vector<double> cov = getCovCpp(rv,x,y,n,rho,tau,eta,p,tol);

  return 0;
}


std::vector<double> getCovCpp(const std::vector<double> & rv,
			      const std::vector<double> & cenX,
			      const std::vector<double> & cenY,
			      const int n,
			      const double rho,
			      const double tau,
			      const double eta,
			      const int p,
			      const double tol){
  Eigen::SparseMatrix<double> sigma;
  std::cout << "sigma" << std::endl;
  getSigmaSparse(sigma,cenX,cenY,n,rho,tau,eta,p,tol);

  std::cout << "decomp" << std::endl;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > llt;
  llt.compute(sigma);

  Eigen::VectorXd covEig;
  Eigen::Map<const Eigen::VectorXd> rvEig(&rv[0],n*p);

  Eigen::SparseMatrix<double> L = llt.matrixL();

  std::cout << "generate" << std::endl;
  covEig = L * rvEig;

  std::vector<double> cov(covEig.data(),
			  covEig.data() + covEig.rows()*covEig.cols());

  std::cout << "done" << std::endl;
  return cov;
}



void getSigmaSparse(Eigen::SparseMatrix<double> & sigma,
		    const std::vector<double> & cenX,
		    const std::vector<double> & cenY,
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
	  dist = std::sqrt(std::pow(cenX.at(i) - cenX.at(k),2.0) +
			   std::pow(cenY.at(i) - cenY.at(k),2.0));
	  val = rho*std::exp(-tau*dist - eta*std::abs(j-l));
	  if(val > tol){
	    sigma.insert(ind0,ind1) = val;
	  }
	}
      }
    }
  }
}
