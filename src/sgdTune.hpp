#ifndef SGD_TUNE_HPP__
#define SGD_TUNE_HPP__




#include <iostream>
#include <omp.h>
#include <map>
#include <algorithm>
#include <limits>
#include "utilities.hpp"
#include "rand.hpp"
#include "data.hpp"
#include "system.hpp"
#include "modelParam.hpp"
#include "model.hpp"
#include "dataDepth.hpp"
#include "agent.hpp"
#include "noTrtAgent.hpp"
#include "myopicAgent.hpp"
#include "proximalAgent.hpp"
#include "rankAgent.hpp"
#include "rankAgentToy.hpp"
#include "runner.hpp"
#include "optim.hpp"
#include "m1SgdOptim.hpp"
#include "m2NmOptim.hpp"



class Design {
 public:

  Design();

  int numVar;
  std::vector<std::string> vars;
  std::map<std::string,int> inds;
  std::map<std::string,double> vals;
  std::map<std::string,std::vector<double> > levs;


  void add(const std::string var, const double val, const double inc,
	   const double lower, const double upper);
  int next();
  void reset();
  void clear();

  double operator [](const std::string) const;
};


#endif
