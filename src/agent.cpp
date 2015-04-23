#include "agent.hpp"


int getNumPre(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  int num =  std::min((int)std::floor(0.06*((double)fD.numNodes)),
		      sD.numNotInfec);
  return std::max(1,num);
}



int getNumAct(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  int num =  std::min((int)std::floor(0.06*((double)fD.numNodes)),
		      sD.numInfected);
  return std::max(1,num);
}
