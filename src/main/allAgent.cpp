#include "allAgent.hpp"

template class AllAgent<ModelGravityGDist>;

template class AllAgent<Model2GravityGDist>;

template class AllAgent<Model2GPowGDist>;

template class AllAgent<ModelGDist>;


template <class M>
std::string AllAgent<M>::name = "all";

template <class M>
void AllAgent<M>::applyTrt(const SimData & sD,
			   TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD,
			   M & m){

  std::for_each(sD.infected.begin(),sD.infected.end(),
		[&tD](const int node){
		  tD.a.at(node) = 1;
		});

  std::for_each(sD.notInfec.begin(),sD.notInfec.end(),
		[&tD](const int node){
		  tD.p.at(node) = 1;
		});

}
