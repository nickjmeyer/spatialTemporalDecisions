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

template<class M, class MP>
class ToyFeatures0 : public BaseFeatures<M,MP> {
 public:
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
  
  std::vector<double> subgraph;

  TrtData tDPre;

  arma::mat infFeat;
  arma::mat notFeat;

  static int numFeatures;

  ToyFeatures0TuneParam tp;

  std::string name;

};


#endif
