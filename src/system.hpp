#ifndef SYSTEM_HPP__
#define SYSTEM_HPP__


#include <vector>
#include <string>
#include "utilities.hpp"
#include "timer.hpp"
#include "data.hpp"
#include "modelParam.hpp"
#include "modelParamGravityTimeInf.hpp"
#include "modelParamGravityTimeInfExp.hpp"
#include "modelParamGravityTimeInfExpCaves.hpp"
#include "modelParamRange.hpp"
#include "modelParamRadius.hpp"
#include "modelParamCave.hpp"
#include "model.hpp"
#include "modelGravityTimeInf.hpp"
#include "modelGravityTimeInfExp.hpp"
#include "modelGravityTimeInfExpCaves.hpp"
#include "modelRange.hpp"
#include "modelRadius.hpp"
#include "modelCave.hpp"
#include "modelMulti.hpp"
#include "rand.hpp"
#include "calcCentrality.hpp"


template <class MG, class ME>
class System {
 public:

  System();
  System(const SimData & sD, const TrtData & tD,
	 const FixedData & fD, const DynamicData & dD,
	 const MG & modelGen, const ME & modelEst);
  System(const std::string file);
  
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;

  SimData sD_r;
  TrtData tD_r;
  DynamicData dD_r;

  MG modelGen;
  ME modelEst;

  MG modelGen_r;
  ME modelEst_r;


  virtual void reset(const std::vector<int> & ind);

  virtual void revert();

  virtual void checkPoint();

  virtual void initialize();

  virtual void preCompData();

  virtual void nextPoint();
  virtual void nextPoint(const std::vector<double> & infProbs);

  virtual void updateStatus();
  
  virtual double value();
};



#endif
