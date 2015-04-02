#ifndef FEATURES_INT_HPP__
#define FEATURES_INT_HPP__

#include "features.hpp"
#include "toyFeatures0.hpp"
#include "toyFeatures1.hpp"
#include "toyFeatures2.hpp"

template <class F, class M>
class FeaturesInt {
 public:
  F f;

  virtual void preCompData(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   M & m);
  
  virtual void getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   M & m);

  virtual void updateFeatures(const SimData & sD,
			      const TrtData & tD,
			      const FixedData & fD,
			      const DynamicData & dD,
			      M & m);

  static int numFeatures;

  arma::mat infFeat;
  arma::mat notFeat;
};


#endif
