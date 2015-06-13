#include "allAgent.hpp"

template class AllAgent<ModelGravityGDist>;

template class AllAgent<ModelTimeExpCavesGDist>;

template class AllAgent<ModelTimeGDistTrendPow>;

template class AllAgent<ModelTimeExpCavesGDistTrendPowCon>;

template class AllAgent<ModelTimeExpCavesEDist>;

template class AllAgent<ModelRadius>;

template class AllAgent<ModelGDist>;

template class AllAgent<ModelCovar>;

template class AllAgent<ModelGDistKern>;


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
