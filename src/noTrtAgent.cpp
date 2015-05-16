#include "noTrtAgent.hpp"


template class NoTrt<ModelGravityGDist>;

template class NoTrt<ModelTimeExpCavesGDist>;

template class NoTrt<ModelTimeGDist>;

template class NoTrt<ModelTimeGDistTrendPow>;

template class NoTrt<ModelTimeExpCavesGDistTrendPowCon>;

template class NoTrt<ModelRadius>;

template<class M>
std::string NoTrt<M>::name = "noTrt";

template<class M>
void NoTrt<M>::applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model){
}

