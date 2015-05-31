#ifndef TUNE_GEN_2_HPP__
#define TUNE_GEN_2_HPP__



#include <iostream>
#include <numeric>
#include <gsl/gsl_multimin.h>
#include "utilities.hpp"
#include "rand.hpp"
#include "data.hpp"
#include "system.hpp"
#include "dataDepth.hpp"
#include "noTrtAgent.hpp"
#include "randomAgent.hpp"
#include "runner.hpp"
#include "settings.hpp"


// typedef ModelTimeExpCavesGDistTrendPowCon M;

template <class M>
class TuneData {
 public:
  TuneData(const System<M,M> s, const std::vector<double> par,
	   const std::vector<double> d, const Starts starts,
	   const int numReps, const double goal0, const double goal1)
    : s(s),par(par),d(d),starts(starts),
      numReps(numReps),goal0(goal0),goal1(goal1){
  }
  System<M,M> s;
  std::vector<double> par;
  std::vector<double> d;
  Starts starts;
  int numReps;
  double goal0;
  double goal1;
};

template <class M>
void scaleD(TuneData<M> * const tuneData, const double & pow);

template <class M>
void tuneSpread(System<M,M> s, const Starts & starts,
		const int numReps,
		double & pow, double & scale,
		const double goal0, const double goal1);

template <class M>
double tuneSpreadHelper(const gsl_vector * x, void * params);

// void tuneTrt(System & s);

// double tuneTrtHelper(const gsl_vector * x, void * params);





#endif
