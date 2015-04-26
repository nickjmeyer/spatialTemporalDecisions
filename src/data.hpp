#ifndef DATA_HPP__
#define DATA_HPP__


#include <vector>
#include <string>
#include <tuple>
#include "settings.hpp"


struct SimData {
  int time; // should also be the length of history

  int numInfected;
  int numNotInfec;

  std::vector<int> infected;
  std::vector<int> notInfec;

  std::vector<int> newInfec;
  std::vector<int> timeInf;
  
  std::vector<int> status;
  std::vector<std::vector<int> > history;

  
  bool operator == (const SimData & rhs) const{
    if(time != rhs.time)
      return false;
    else if(numInfected != rhs.numInfected)
      return false;
    else if(numNotInfec != rhs.numNotInfec)
      return false;
    else if(infected != rhs.infected)
      return false;
    else if(notInfec != rhs.notInfec)
      return false;
    else if(newInfec != rhs.newInfec)
      return false;
    else if(timeInf != rhs.timeInf)
      return false;
    else if(status != rhs.status)
      return false;
    else if(history != rhs.history)
      return false;
    else
      return true;
  };
};


struct TrtData {
  std::vector<int> a;
  std::vector<int> p;

  std::vector<int> aPast;
  std::vector<int> pPast;


  bool operator == (const TrtData & rhs) const {
    if(a != rhs.a)
      return false;
    else if(p != rhs.p)
      return false;
    else if(aPast != rhs.aPast)
      return false;
    else if(pPast != rhs.pPast)
      return false;
    else
      return true;
  };
};


struct FixedData {
  int numNodes;

  int trtStart;

  int period;

  int finalT;

  // double propTrt;

  double priorTrtMean;
  
  std::vector<double> dist;

  std::vector<double> caves;

  std::vector<double> covar;
  int numCovar;

  std::vector<int> fips;

  std::vector<int> network;
  
  std::vector<double> centroidsLong;
  std::vector<double> centroidsLat;

  std::vector<double> centroidsMdsLong;
  std::vector<double> centroidsMdsLat;
  
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

  double distSD;
  std::vector<double> expDistSD; // e^{-dist^2/(2*sd(dist)^2)}
  
  std::vector<double> logDist; // log(2+dist)

  std::vector<double> hpdd;
};


struct DynamicData {

  bool operator == (const DynamicData & rhs) const {
    return true;
  };
};


typedef std::tuple<SimData,TrtData,DynamicData> DataBundle;
std::vector<DataBundle>
historyToData(const std::vector<std::vector<int> > & hist);



#endif
