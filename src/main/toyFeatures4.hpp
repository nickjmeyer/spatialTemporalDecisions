#ifndef TOY_FEATURES_4_HPP__
#define TOY_FEATURES_4_HPP__

#include <vector>
#include <queue>
#include "features.hpp"
#include "tuneParam.hpp"
#include "calcCentrality.hpp"

class ToyFeatures4TuneParam : public TuneParam {
 public:
  virtual std::vector<double> getPar() const ;
  virtual void putPar(const std::vector<double> & par);
};

template<class M>
class ToyFeatures4 : public BaseFeatures<M> {
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
  std::vector<std::vector<std::pair<int,double> > > notNeigh;
  std::vector<int> notNeighNum;

  // not infected are neighbors of
  std::vector<std::vector<std::pair<int,int> > > notNeighOf; 
  std::vector<int> notNeighOfNum;

  // sub graph of not infec
  arma::colvec subGraphNotInfec;

  TrtData tDPre;

  arma::mat infFeat;
  arma::mat notFeat;

  const static int numFeatures = 4;

  ToyFeatures4TuneParam tp;

  std::string name;

};



#endif
