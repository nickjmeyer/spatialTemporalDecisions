#ifndef TUNE_SP_HPP__
#define TUNE_SP_HPP__


#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include "utilities.hpp"
#include "rand.hpp"
#include "data.hpp"
#include "system.hpp"
#include "modelParam.hpp"
#include "model.hpp"
#include "modelParamEbola.hpp"
#include "modelEbola.hpp"
#include "dataDepth.hpp"
#include "noTrtAgent.hpp"
#include "myopicAgent.hpp"
#include "rankAgentToy.hpp"
#include "features.hpp"
#include "toyFeatures2.hpp"
#include "featuresInt.hpp"
#include "runner.hpp"
#include "settings.hpp"
#include "m1SpOptim.hpp"



class FFX {
 public:

  FFX();

  std::vector<std::string> factors;
  std::vector<std::vector<double> > values;
  
  std::vector<std::string> stats;

  std::vector<std::vector<double> > results;

  std::vector<int> maxSett;
  
  void addFactor(const std::string & f,
		 const std::vector<double> & fVals);
  void addStat(const std::string & s);

  int numFactor;

  int numCombo;
  int numReps;

  int getMax() const ;
  std::vector<double> getSett(const int i) const;

  
  
};



#endif
