#ifndef FEATURES_HPP__
#define FEATURES_HPP__

#include "data.hpp"
#include "system.hpp"

template <class M, class MP>
class BaseFeatures {
 public:
  virtual void preCompData(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const M & m,
			   MP & mP) = 0;
  
  virtual void getFeatures(const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const M & m,
			   MP & mP) = 0;

  virtual void updateFeatures(const SimData & sD,
			      const TrtData & tD,
			      const FixedData & fD,
			      const DynamicData & dD,
			      const M & m,
			      MP & mP) = 0;

  TrtData tDPre;

  arma::mat infFeat;
  arma::mat notFeat;

};




void reconstructNetwork(std::vector<int> & network,
			std::vector<int> & treated);




#endif
