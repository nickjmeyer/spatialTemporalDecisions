#include "agent.hpp"


int getNumPre(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  return std::min((int)std::floor(0.03*((double)fD.numNodes)+1),
  		  sD.numNotInfec);
  // return std::min(1,sD.numNotInfec);
}



int getNumAct(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  return std::min((int)std::floor(0.03*((double)fD.numNodes)+1),
  		  sD.numInfected);
  // return std::min(1,sD.numInfected);
}
