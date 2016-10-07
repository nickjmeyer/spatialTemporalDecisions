#include "agent.hpp"


int getNumPre(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  // int num =  std::min((int)std::floor(0.06*((double)fD.numNodes)),
  // 		      sD.numNotInfec);
  // return std::max(1,num);
  int num =  std::max((int)std::floor(0.06*((double)fD.numNodes)),
  		      1);
  return std::min(sD.numNotInfec,num);
}



int getNumAct(const SimData & sD,
	      const TrtData & tD,
	      const FixedData & fD,
	      const DynamicData & dD){
  // int num =  std::min((int)std::floor(0.06*((double)fD.numNodes)),
  // 		      sD.numInfected);
  // return std::max(1,num);
  int num = std::max((int)std::floor(0.06*((double)fD.numNodes)),
  		     1);
  return std::min(sD.numInfected,num);
}


void checkForValidTrt(
  const SimData & sD, const TrtData & tD,
    const FixedData & fD, const DynamicData & dD) {
  const int numAct = getNumAct(sD,tD,fD,dD);
  const int numPre = getNumPre(sD,tD,fD,dD);

  int totPre = 0,totAct = 0;
  // check if valid treatments are given to valid locations
  for(int i = 0; i < fD.numNodes; i++){
    CHECK(tD.p.at(i) == 1 || tD.p.at(i) == 0)
      << "Prevenative treatment not 1 or 0"
      << ": " << tD.p.at(i)
      << std::endl;
    CHECK(tD.a.at(i) == 1 || tD.a.at(i) == 0)
      << "Active treatment not 1 or 0"
      << std::endl;
    if(tD.a.at(i) == 1) {
      CHECK(sD.status.at(i) >= 2)
        << "Not infected receiving active treatment"
        << std::endl;
    }
    if(tD.p.at(i) == 1) {
      CHECK(sD.status.at(i) >= 2)
        << "Infected receiving preventative treament"
        << std::endl;
    }
    else if(tD.a.at(i) == 1)
      totAct++;
    else if(tD.p.at(i) == 1)
      totPre++;
  }

  // check if total number of treatments are correct
  CHECK(totAct == numAct)
    << "Not correct amount of active treatments."
    << std::endl
    << "Should be " << numAct << " but is " << totAct << "."
    << std::endl
    << "Number of infected nodes is " << sD.numInfected
    << "(" << sD.infected.size() << ")"
    << std::endl;
  CHECK(totPre == numPre)
    << "Not correct amount of preventative treatments."
    << std::endl
    << "Should be " << numPre << " but is " << totPre << "."
    << std::endl
    << "Number of not infected nodes is " << sD.numNotInfec
    << "(" << sD.notInfec.size() << ")"
    << std::endl;
}
