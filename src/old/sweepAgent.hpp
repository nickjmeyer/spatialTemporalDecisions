#ifndef SWEEP_AGENT_HPP__
#define SWEEP_AGENT_HPP__


#include <armadillo>
#include <vector>
#include <queue>
#include <limits>
#include <numeric>
#include "system.hpp"
#include "agent.hpp"
#include "dataDepth.hpp"
#include "tuneParam.hpp"



class SweepTuneParam : public TuneParam {
 public:
  std::vector<double> getPar() const;
  void putPar(const std::vector<double> & par);

  double sigma;
  arma::colvec weights;

  int numChunks;
  int maxSweep;
};




template<class Model, class ModelParam>
class SweepAgent : public BaseAgent<Model,ModelParam>{
 public:
  SweepAgent();
  
  void applyTrt(const SimData & sD,
		TrtData & tD,
		const FixedData & fD,
		const DynamicData & dD,
		const Model & m,
		ModelParam & mP);

  virtual void getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   const ModelParam & mP);

  virtual double calcValue();

  virtual void updateSameP(const int ind, const int change,
			   const SimData & sD,
			   TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);
  virtual void updateSwapP(const int oldInd, const int newInd,
			   const SimData & sD,
			   TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);
  virtual void updateSameA(const int ind, const int change,
			   const SimData & sD,
			   TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);
  virtual void updateSwapA(const int oldInd, const int newInd,
			   const SimData & sD,
			   TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);
  
  
  arma::mat infFeat;
  arma::mat notFeat;
  arma::colvec infRanks;
  arma::colvec notRanks;
  
  static int numFeatures;

  int numAct;
  int numPre;

  SweepTuneParam tp;

  std::string name;

};


#endif
