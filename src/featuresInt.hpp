#ifndef FEATURES_INT_HPP__
#define FEATURES_INT_HPP__

#include "features.hpp"
#include "toyFeatures0.hpp"
#include "toyFeatures1.hpp"
#include "toyFeatures2.hpp"

template <class F, class M, class MP>
class FeaturesInt {
 public:
  F f;

  virtual void preCompData(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const M & m,
			   MP & mP);
  
  virtual void getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const M & m,
			   MP & mP);

  virtual void updateFeatures(const SimData & sD,
			      const TrtData & tD,
			      const FixedData & fD,
			      const DynamicData & dD,
			      const M & m,
			      MP & mP);

  static int numFeatures;

  arma::mat infFeat;
  arma::mat notFeat;
};


#endif
