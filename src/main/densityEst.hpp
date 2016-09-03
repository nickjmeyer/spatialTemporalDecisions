#ifndef DENSITY_EST_HPP__
#define DENSITY_EST_HPP__

#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <queue>
#include <gsl/gsl_multimin.h>


class DensityEst {
 private:
  double gKernel(const double x);
  double gKernelP(const double x);

  static double f(const gsl_vector * x, void * params);
  static void df(const gsl_vector * x, void * params,
		 gsl_vector * g);
  static void fdf(const gsl_vector * x, void * params,
		  double * f, gsl_vector * g);

  double optH(const std::vector<double> x);

  std::vector<double> xi;
  double h;

 public:
  DensityEst(const std::vector<double> xi);

  double eval(const double x);
  double deriv(const double x);

  std::pair<double,double> max();
};




#endif
