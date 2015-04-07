#include "agent.hpp"


int getNumPre(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  int num =  std::min((int)std::floor(fD.propTrt*((double)fD.numNodes)+1),
		      sD.numNotInfec);
  return std::max(1,num);
}



int getNumAct(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  int num =  std::min((int)std::floor(fD.propTrt*((double)fD.numNodes)+1),
		      sD.numInfected);
  return std::max(1,num);
}
