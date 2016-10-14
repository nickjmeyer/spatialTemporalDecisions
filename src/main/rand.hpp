#ifndef RAND_HPP
#define RAND_HPP

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <vector>

namespace njm{
  void resetSeed();
  void resetSeed(const int seed);
  void resetSeedAll();


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


typedef boost::variate_generator<boost::mt19937,
				 boost::uniform_real<> > VGRunif01;

typedef boost::variate_generator<boost::mt19937,
				 boost::normal_distribution<> > VGRnorm01;

class RandParr{
 public:

  RandParr(const int numSource, const int numRand);

 protected:

  void fillRunif01(const int source);
  void fillRnorm01(const int source);

  void reset();
  void reset(const int source);

  void setSeed(const int source, const int seed);

  double genRunif01(const int source);
  double genRnorm01(const int source);


  int numSource;
  int numRand;

  std::vector<int> seeds;
  std::vector<VGRunif01> vgRunif01;
  std::vector<VGRnorm01> vgRnorm01;

  std::vector< std::vector<double>::iterator > runif01Iter;
  std::vector< std::vector<double>::iterator > runif01End;
  std::vector< std::vector<double> > runif01Vals;

  std::vector< std::vector<double>::iterator > rnorm01Iter;
  std::vector< std::vector<double>::iterator > rnorm01End;
  std::vector< std::vector<double> > rnorm01Vals;


  friend void njm::resetSeed();
  friend void njm::resetSeed(const int seed);
  friend void njm::resetSeedAll();

  friend double njm::runif01();
  friend double njm::rnorm01();
};



#endif
