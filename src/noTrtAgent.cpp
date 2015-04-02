#include "noTrtAgent.hpp"


template class NoTrt<GravityModel>;
template class NoTrt<GravityTimeInfModel>;
template class NoTrt<GravityTimeInfExpCavesModel>;

template class NoTrt<RangeModel>;

template class NoTrt<RadiusModel>;

template class NoTrt<CaveModel>;

template<class M>
std::string NoTrt<M>::name = "noTrt";

template<class M>
void NoTrt<M>::applyTrt(const SimData & sD,
			TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			M & model){
}

