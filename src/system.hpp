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


template <class Model, class ModelParam>
class System {
 public:

  System();
  System(const SimData & sD, const TrtData & tD,
	 const FixedData & fD, const DynamicData & dD,
	 const Model & model,
	 const ModelParam & genParam, const ModelParam & estParam);
  
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;

  SimData sD_r;
  TrtData tD_r;
  DynamicData dD_r;

  Model model;
  ModelParam genParam;
  ModelParam estParam;

  ModelParam genParam_r;
  ModelParam estParam_r;


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
	      const ModelParam & genParam);
  
  SimData sD;
  TrtData tD;
  FixedData fD;
  DynamicData dD;

  SimData sD_r;
  TrtData tD_r;
  DynamicData dD_r;

  Model model;
  ModelParam genParam;

  ModelParam genParam_r;


  virtual void reset();

  virtual void nextPoint(const int isFinal = 0);
  
  virtual double value();
};


#endif
