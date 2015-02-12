#include "noTrtAgent.hpp"


template class NoTrt<GravityModel,GravityParam>;

template class NoTrt<RangeModel,RangeParam>;

template class NoTrt<EbolaModel,EbolaParam>;

template<class M, class MP>
std::string NoTrt<M,MP>::name = "noTrt";

template<class M, class MP>
void NoTrt<M,MP>::applyTrt(const SimData & sD,
			   TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   const M & model,
			   MP & modelParam){
}

