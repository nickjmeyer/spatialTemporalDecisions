#include "noTrtAgent.hpp"


template class NoTrt<ModelGravity>;

template class NoTrt<ModelTimeExpCaves>;

template class NoTrt<ModelTime>;

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

