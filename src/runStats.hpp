#ifndef RUN_STATS_HPP__
#define RUN_STATS_HPP__

#include <vector>
#include <cmath>
#include <algorithm>

class RunStats {
 protected:
  std::vector<double> vals;
  unsigned int n;

  double smean_; // sample mean
  double svar_; // sample variance
  double ssd_; // sample standard deviation
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

  double smean() const;

  double svar() const;

  double ssd() const;

  double seMean() const;
};




#endif
