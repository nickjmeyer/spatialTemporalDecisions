#ifndef SYSTEM_HPP__
#define SYSTEM_HPP__


#include <vector>
#include <string>
#include "utilities.hpp"
#include "data.hpp"
#include "modelParam.hpp"
#include "modelParamGravityTimeInf.hpp"
#include "modelParamGravityTimeInfSq.hpp"
#include "modelParamGravityTimeInfSqrt.hpp"
#include "modelParamGravityTimeInfLog.hpp"
#include "modelParamGravityTimeInfExp.hpp"
#include "modelParamGravityTimeInfExpCaves.hpp"
#include "modelParamGravityTimeInfExpLCaves.hpp"
#include "modelParamGravityTimeInfExpRCaves.hpp"
#include "modelParamEbola.hpp"
#include "modelParamRange.hpp"
#include "modelParamRadius.hpp"
#include "modelParamCave.hpp"
#include "model.hpp"
#include "modelGravityTimeInf.hpp"
#include "modelGravityTimeInfSq.hpp"
#include "modelGravityTimeInfSqrt.hpp"
#include "modelGravityTimeInfLog.hpp"
#include "modelGravityTimeInfExp.hpp"
#include "modelGravityTimeInfExpCaves.hpp"
#include "modelGravityTimeInfExpLCaves.hpp"
#include "modelGravityTimeInfExpRCaves.hpp"
#include "modelEbola.hpp"
#include "modelRange.hpp"
#include "modelRadius.hpp"
#include "modelCave.hpp"
#include "rand.hpp"
#include "calcCentrality.hpp"


template <class MG, class MPG,
	  class ME, class MPE>
class System {
 public:

  System();
  System(const SimData & sD, const TrtData & tD,
	 const FixedData & fD, const DynamicData & dD,
	 const MG & modelGen, const ME & modelEst,
	 const MPG & paramGen, const MPE & paramEst);
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
  MPG paramGen;
  MPE paramEst;

  MPG paramGen_r;
  MPE paramEst_r;


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



template <class M, class MP>
class SystemLight {
 public:

  SystemLight(const SimData & sD, const TrtData & tD,
	      const FixedData & fD, const DynamicData & dD,
	      const M & model,
	      const MP & paramGen);
  
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;

  SimData sD_r;
  TrtData tD_r;
  DynamicData dD_r;

  M modelGen;
  MP paramGen;

  MP paramGen_r;


  virtual void reset();

  virtual void nextPoint(const int isFinal = 0);
  
  virtual double value();
};


#endif
