#ifndef TOY_FEATURES_2_MULTI_HPP__
#define TOY_FEATURES_2_MULTI_HPP__

#include <vector>
#include <queue>
#include "features.hpp"
#include "tuneParam.hpp"
#include "calcCentrality.hpp"
#include "modelMulti.hpp"

class ToyFeatures2MultiTuneParam : public TuneParam {
 public:
  virtual std::vector<double> getPar() const ;
  virtual void putPar(const std::vector<double> & par);
};

template<class M>
class ToyFeatures2Multi : public BaseFeatures<M> {
 public:
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
  
  // neighbors of not infected
  // models[ uninfected[ neighbors and probs[...] ] ]
  std::vector<std::vector<std::vector<std::pair<int,double> > > > notNeigh;
  std::vector<int> notNeighNum;

  // not infected are neighbors of
  std::vector<std::vector<std::pair<int,int> > > notNeighOf; 
  std::vector<int> notNeighOfNum;

  // sub graph of not infec
  arma::colvec subGraphNotInfec;

  TrtData tDPre;

  arma::mat infFeat;
  arma::mat notFeat;

  const static int numFeatures = M::numModels*3 + 1;

  ToyFeatures2MultiTuneParam tp;

  std::string name;

};



#endif
