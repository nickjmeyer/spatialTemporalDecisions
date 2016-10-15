#include "randomAgent.hpp"


template<class M>
const std::string RandomAgent<M>::name = "random";

template <class M>
void RandomAgent<M>::applyTrt(const SimData & sD,
        TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD,
        M & m){
    numPre = getNumPre(sD,tD,fD,dD);
    numAct = getNumAct(sD,tD,fD,dD);


    std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;

    int i,node0;
    for(i=0; i<sD.numNotInfec; i++){
        node0=sD.notInfec.at(i);
        sortNotInfec.push(std::pair<double,int>(njm::runif01(),node0));
    }


    for(i=0; i<sD.numInfected; i++){
        node0=sD.infected.at(i);
        sortInfected.push(std::pair<double,int>(njm::runif01(),node0));
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

template class RandomAgent<ModelGravityGDist>;
