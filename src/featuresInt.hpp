#ifndef FEATURES_INT_HPP__
#define FEATURES_INT_HPP__

#include "features.hpp"
#include "toyFeatures0.hpp"
#include "toyFeatures1.hpp"
#include "toyFeatures2.hpp"

template <class Features, class Model, class ModelParam>
class FeaturesInt {
 public:
  Features f;

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

  static int numFeatures;

  arma::mat infFeat;
  arma::mat notFeat;
};


#endif
