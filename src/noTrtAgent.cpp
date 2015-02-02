#include "noTrtAgent.hpp"


template class NoTrt<GravityModel,GravityParam>;

template class NoTrt<RangeModel,RangeParam>;

template class NoTrt<EbolaModel,EbolaParam>;

template<class Model, class ModelParam>
const std::string NoTrt<Model,ModelParam>::name = "noTrt";

template<class Model, class ModelParam>
void NoTrt<Model,ModelParam>::applyTrt(const SimData & sD,
				       TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD,
				       const Model & model,
				       ModelParam & modelParam){
}

