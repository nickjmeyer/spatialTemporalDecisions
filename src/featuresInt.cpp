#include "featuresInt.hpp"

template class FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
			   GravityModel,GravityParam>;

template class FeaturesInt<ToyFeatures2<RangeModel,RangeParam>,
			   RangeModel,RangeParam>;

template class FeaturesInt<ToyFeatures2<CaveModel,CaveParam>,
			   CaveModel,CaveParam>;


template <class F, class M, class MP>
int FeaturesInt<F,M,MP>::numFeatures = F::numFeatures + 1;

template <class F, class M, class MP>
void FeaturesInt<F,M,MP>::preCompData(const SimData & sD,
				      const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD,
				      const M & m,
				      MP & mP){
  f.preCompData(sD,tD,fD,dD,m,mP);
}


template <class F, class M, class MP>
void FeaturesInt<F,M,MP>::getFeatures(const SimData & sD,
				      const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD,
				      const M & m,
				      MP & mP){
  f.getFeatures(sD,tD,fD,dD,m,mP);
  
  infFeat.ones(sD.numInfected,numFeatures);
  notFeat.ones(sD.numNotInfec,numFeatures);

  infFeat.submat(0,1,infFeat.n_rows-1,infFeat.n_cols-1) = f.infFeat;
  notFeat.submat(0,1,notFeat.n_rows-1,notFeat.n_cols-1) = f.notFeat;
}


template <class F, class M, class MP>
void
FeaturesInt<F,M,MP>::updateFeatures(const SimData & sD,
				    const TrtData & tD,
				    const FixedData & fD,
				    const DynamicData & dD,
				    const M & m,
				    MP & mP){
  f.updateFeatures(sD,tD,fD,dD,m,mP);
  
  infFeat.ones(sD.numInfected,numFeatures);
  notFeat.ones(sD.numNotInfec,numFeatures);

  infFeat.submat(0,1,infFeat.n_rows-1,infFeat.n_cols-1) = f.infFeat;
  notFeat.submat(0,1,notFeat.n_rows-1,notFeat.n_cols-1) = f.notFeat;
}

