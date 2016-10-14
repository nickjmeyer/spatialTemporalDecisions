#ifndef TUNE_TRT_PROP_HPP
#define TUNE_TRT_PROP_HPP


#include <iostream>
#include <vector>
#include <string>
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
#include "rankAgent.hpp"
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

  std::vector<std::vector<double> > allObs;

  std::vector<int> maxSett;

  void setReps(const int num);

  void addFactor(const std::string & f,
		 const std::vector<double> & fVals);
  void addStat(const std::string & s);

  int numFactor;
  int numStat;

  int numCombo;
  int numReps;

  int maxInd;

  int getMax() const ;
  std::vector<double> getSett(const int num) const;
  double getSett(const std::string & f, const int num) const;

  void addObs(const int num, const double & obs);
  void addObs(const int num, const std::vector<double> & obs);

  void saveObs(const std::string & file) const;
};



#endif
