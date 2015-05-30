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


typedef ModelTimeExpCavesGDistTrendPowCon M;

class TuneData {
 public:
  TuneData(System<M,M> s, std::vector<double> par,
	   std::vector<double> d, Starts starts, int numReps,
	   double goal)
    : s(s),par(par),d(d),starts(starts),numReps(numReps),goal(goal){
  }
  System<M,M> s;
  std::vector<double> par;
  std::vector<double> d;
  Starts starts;
  int numReps;
  double goal;
};

void scaleD(TuneData<M> * const tuneData, const double & pow);

void tuneSpread(System<M,M> s, const Starts & starts,
		const int numReps, double & pow,
		double & scale, const double goal);

double tuneSpreadHelper(const gsl_vector * x, void * params);

// void tuneTrt(System & s);

// double tuneTrtHelper(const gsl_vector * x, void * params);





#endif
