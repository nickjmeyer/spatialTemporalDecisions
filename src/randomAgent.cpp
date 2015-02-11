#include "randomAgent.hpp"


template class RandomAgent<GravityModel,GravityParam>;

template<class M, class MP>
const std::string RandomAgent<M,MP>::name = "random";

template <class M, class MP>
void RandomAgent<M,MP>::applyTrt(const SimData & sD,
				 TrtData & tD,
				 const FixedData & fD,
				 const DynamicData & dD,
				 const M & m,
				 MP & mP){
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);


  std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;

  for(i=0; i<sD.numNotInfec; i++){
    node0=sD.notInfec.at(i);
    sortNotInfec.push(std::pair<double,int>(njm::runif(),node0));
  }    

  
  for(i=0; i<sD.numInfected; i++){
    node0=sD.infected.at(i);
    sortInfected.push(std::pair<double,int>(njm::runif(),node0));
  }    
	 

  for(i=0; i<numAct; i++){
    tD.a.at(sortInfected.top().second) = 1;
    sortInfected.pop();
  }

  for(i=0; i<numPre; i++){
    tD.p.at(sortNotInfec.top().second) = 1;
    sortNotInfec.pop();
  }
}

