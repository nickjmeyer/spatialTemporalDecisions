#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
// #include <eigen3/Eigen/Eigen>
// #include <eigen3/Eigen/Sparse>
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
arma::mat getSigmaSparse(const std::vector<double> & cenX,
        const std::vector<double> & cenY,
        const int n,
        const double rho,
        const double tau,
        const double eta,
        const int p);


std::vector<double> getCovCpp(const std::vector<double> & rv,
        const std::vector<double> & cenX,
        const std::vector<double> & cenY,
        const int n,
        const double rho,
        const double tau,
        const double eta,
        const int p){
    arma::mat sigma;
    sigma = getSigmaSparse(cenX,cenY,n,rho,tau,eta,p);


    arma::colvec rvArma(rv);

    arma::mat r = arma::chol(sigma);

    return arma::conv_to<std::vector<double> >::from(r.t()*rvArma);
}



arma::mat getSigmaSparse(const std::vector<double> & cenX,
        const std::vector<double> & cenY,
        const int n,
        const double rho,
        const double tau,
        const double eta,
        const int p){
    int i,j,k,l,ind0,ind1;
    double val;
    double dist;

    arma::mat sigma;
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
    return(sigma);
}
