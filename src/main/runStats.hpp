#ifndef RUN_STATS_HPP
#define RUN_STATS_HPP

#include <vector>
#include <cmath>
#include <algorithm>

class RunStats {
 protected:
  std::vector<double> vals;
  unsigned int n;

  double sMean_; // sample mean
  double sSqMean_; // sample mean of squared values
  double sVar_; // sample variance
  double sSd_; // sample standard deviation
  double seMean_; // standard error of mean

  void update(const double & add);
  void update(const std::vector<double> & add);

 public:
  RunStats();
  RunStats(const std::vector<double> & init);

  void operator () (const std::vector<double> & add){
    update(add);
  };

  void operator () (const double & add){
    update(add);
  };

  double sMean() const;

  double sVar() const;

  double sSd() const;

  double seMean() const;
};




#endif
