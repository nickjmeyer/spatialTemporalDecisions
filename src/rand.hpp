#ifndef RAND_HPP__
#define RAND_HPP__

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <vector>

class RandParr{
 public:

  RandParr(int numRand_);
  void reset();
  void fixSeed(const int fix);
  double getRunif01();
  double getRnorm01();


  int numRand;

 private:

  void initialize();

  int numRunif01;
  std::vector< std::vector<double>::iterator > runif01Iter, runif01Iter_hold;
  std::vector< std::vector<double>::iterator > runif01End, runif01End_hold;
  std::vector< std::vector<double> > runif01Vals, runif01Vals_fixed;
  
  int numRnorm01;
  std::vector< std::vector<double>::iterator > rnorm01Iter, rnorm01Iter_hold;
  std::vector< std::vector<double>::iterator > rnorm01End, rnorm01End_hold;
  std::vector< std::vector<double> > rnorm01Vals, rnorm01Vals_fixed;

  std::vector< bool > isfixed;
};


void resetRandomSeed();
void fixRandomSeed(const int fix);

namespace njm{
  double runif01();
  double runif(double a, double b);
  
  int runifInterv(int min, int max);
  // REQUIRES: max, min are positive integers, max>min
  // MODIFIES: nothing
  // EFFECTS: returns a random integer from the set {0,1,...,max-1}
  
  int rber(double p);
  // REQUIRES: 0<=p<=1
  // MODIFIES: nothing
  // EFFECTS: returns 1 with probability p and 0 with probability (1-p)   
  
  double rnorm01();
  // REQUIRES: nothing
  // MODIFIES: nothing
  // EFFECTS: returns a sample from a standard normal distribution
  
  double rgamma(double const alpha, double const beta);
  
  double qnormal(double prob);
  double pnormal(double quan);
};

#endif
