#ifndef DATA_HPP__
#define DATA_HPP__


#include <vector>
#include <string>
#include "settings.hpp"


struct SimData {
  int time;

  int numInfected;
  int numNotInfec;

  std::vector<int> infected;
  std::vector<int> notInfec;

  std::vector<int> newInfec;
  std::vector<int> timeInf;
  
  std::vector<int> status;
  std::vector<std::vector<int> > history;
};


struct TrtData {
  std::vector<int> a;
  std::vector<int> p;

  std::vector<int> aPast;
  std::vector<int> pPast;
};


struct FixedData {
  int numNodes;

  int trtStart;

  int period;

  int finalT;

  double priorTrtMean;
  
  std::vector<double> dist;

  std::vector<double> caves;

  std::vector<double> covar;
  int numCovar;

  std::vector<int> fips;

  std::vector<int> network;
  
  std::vector<double> centroidsLong;
  std::vector<double> centroidsLat;

  std::vector<double> subGraph;
  std::vector<double> betweenness;


  
  // pre-computed data
  std::vector<double> propCaves;
  std::vector<double> logPropCaves;
  std::vector<double> rankCaves;
  
  std::vector<double> subGraphK;
  int subGraphKval;
  double subGraphKmax;
  
  double invDistSD;
  std::vector<double> expInvDistSD; // e^{[1/(1+dist)]/[sd(1+dist)]}
  std::vector<double> logDist; // log(2+dist)
  
};


struct DynamicData {
};


#endif
