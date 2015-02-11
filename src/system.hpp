#ifndef SYSTEM_HPP__
#define SYSTEM_HPP__


#include <vector>
#include <string>
#include "utilities.hpp"
#include "data.hpp"
#include "modelParam.hpp"
#include "modelParamEbola.hpp"
#include "modelParamRange.hpp"
#include "model.hpp"
#include "modelEbola.hpp"
#include "modelRange.hpp"
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


  virtual void reset();

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
