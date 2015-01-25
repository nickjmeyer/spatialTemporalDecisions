#include "featuresInt.hpp"

template class FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
			   GravityModel,GravityParam>;

template <class Features, class Model, class ModelParam>
int FeaturesInt<Features,Model,ModelParam>::numFeatures =
  Features::numFeatures + 1;

template <class Features, class Model, class ModelParam>
void FeaturesInt<Features,Model,ModelParam>::preCompData(const SimData & sD,
							 const TrtData & tD,
							 const FixedData & fD,
							 const DynamicData & dD,
							 const Model & m,
							 ModelParam & mP){
  f.preCompData(sD,tD,fD,dD,m,mP);
}


template <class Features, class Model, class ModelParam>
void FeaturesInt<Features,Model,ModelParam>::getFeatures(const SimData & sD,
							 const TrtData & tD,
							 const FixedData & fD,
							 const DynamicData & dD,
							 const Model & m,
							 ModelParam & mP){
  f.getFeatures(sD,tD,fD,dD,m,mP);
  
  infFeat.ones(sD.numInfected,numFeatures);
  notFeat.ones(sD.numNotInfec,numFeatures);

  infFeat.submat(0,1,infFeat.n_rows-1,infFeat.n_cols-1) = f.infFeat;
  notFeat.submat(0,1,notFeat.n_rows-1,notFeat.n_cols-1) = f.notFeat;
}


template <class Features, class Model, class ModelParam>
void
FeaturesInt<Features,Model,ModelParam>::updateFeatures(const SimData & sD,
						       const TrtData & tD,
						       const FixedData & fD,
						       const DynamicData & dD,
						       const Model & m,
						       ModelParam & mP){
  f.updateFeatures(sD,tD,fD,dD,m,mP);
  
  infFeat.ones(sD.numInfected,numFeatures);
  notFeat.ones(sD.numNotInfec,numFeatures);

  infFeat.submat(0,1,infFeat.n_rows-1,infFeat.n_cols-1) = f.infFeat;
  notFeat.submat(0,1,notFeat.n_rows-1,notFeat.n_cols-1) = f.notFeat;
}

