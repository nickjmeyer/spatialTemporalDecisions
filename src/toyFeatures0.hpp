#ifndef TOY_FEATURES_0_HPP__
#define TOY_FEATURES_0_HPP__

#include <vector>
#include <queue>
#include "features.hpp"
#include "tuneParam.hpp"
#include "calcCentrality.hpp"

class ToyFeatures0TuneParam : public TuneParam {
 public:
  virtual std::vector<double> getPar() const ;
  virtual void putPar(const std::vector<double> & par);

  int valReps;
};

template<class Model, class ModelParam>
class ToyFeatures0 : public BaseFeatures<Model,ModelParam> {
 public:
  virtual void preCompData(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);

  virtual void getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const Model & m,
			   ModelParam & mP);

  virtual void updateFeatures(const SimData & sD,
			      const TrtData & tD,
			      const FixedData & fD,
			      const DynamicData & dD,
			      const Model & m,
			      ModelParam & mP);
  
  std::vector<double> subgraph;

  TrtData tDPre;

  arma::mat infFeat;
  arma::mat notFeat;

  static int numFeatures;

  ToyFeatures0TuneParam tp;

  std::string name;

};


#endif
