#include "noTrtAgent.hpp"


template<class M>
std::string NoTrt<M>::name = "noTrt";

template<class M>
void NoTrt<M>::applyTrt(const SimData & sD,
        TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD,
        M & model){
}


template class NoTrt<ModelGravityGDist>;

template class NoTrt<Model2GravityGDist>;

template class NoTrt<Model2GravityEDist>;

template class NoTrt<Model2GPowGDist>;

template class NoTrt<Model2EdgeToEdge>;

template class NoTrt<ModelGDist>;
