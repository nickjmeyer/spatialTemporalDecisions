#ifndef SYSTEM_HPP__
#define SYSTEM_HPP__


#include <vector>
#include <string>
#include "utilities.hpp"
#include "data.hpp"
#include "modelParam.hpp"
#include "model.hpp"
#include "modelParamEbola.hpp"
#include "modelEbola.hpp"
#include "rand.hpp"
#include "calcCentrality.hpp"


template <class ModelGen, class ModelParamGen,
	  class ModelEst, class ModelParamEst>
class System {
 public:

  System();
  System(const SimData & sD, const TrtData & tD,
	 const FixedData & fD, const DynamicData & dD,
	 const ModelGen & modelGen, const ModelEst & modelEst,
	 const ModelParamGen & paramGen, const ModelParamEst & paramEst);
  
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;

  SimData sD_r;
  TrtData tD_r;
  DynamicData dD_r;

  ModelGen modelGen;
  ModelEst modelEst;
  ModelParamGen paramGen;
  ModelParamEst paramEst;

  ModelParamGen paramGen_r;
  ModelParamEst paramEst_r;


  virtual void reset();

  virtual void checkPoint();

  virtual void initialize();

  virtual void preCompData();

  virtual void nextPoint();
  virtual void nextPoint(const std::vector<double> & infProbs);

  virtual void updateStatus();
  
  virtual double value();
};



template <class Model, class ModelParam>
class SystemLight {
 public:

  SystemLight(const SimData & sD, const TrtData & tD,
	      const FixedData & fD, const DynamicData & dD,
	      const Model & model,
	      const ModelParam & paramGen);
  
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;

  SimData sD_r;
  TrtData tD_r;
  DynamicData dD_r;

  Model modelGen;
  ModelParam paramGen;

  ModelParam paramGen_r;


  virtual void reset();

  virtual void nextPoint(const int isFinal = 0);
  
  virtual double value();
};


#endif
