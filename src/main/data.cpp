#include "data.hpp"


std::vector<DataBundle>
historyToData(const std::vector<std::vector<int> > & hist){
  int years = int(hist.size()), numNodes = int(hist[0].size());
  std::vector<DataBundle> res;
  
  SimData sD;
  TrtData tD;
  DynamicData dD;
  std::vector<int> aPast,pPast;

  sD.timeInf.resize(numNodes);
  std::fill(sD.timeInf.begin(),sD.timeInf.end(),0);
  
  int i,j,hij;
  for(i = 0; i < years; ++i){
    sD.time = i;

    // wipe everything
    sD.numInfected = 0;
    sD.numNotInfec = 0;
    sD.infected.clear();
    sD.notInfec.clear();
    sD.newInfec.clear();
    // update history before wiping
    if(i > 0)
      sD.history.push_back(sD.status);
    sD.status.clear();


    // grab past treatment before wiping
    if(i == 0){
      tD.aPast.resize(numNodes);
      tD.pPast.resize(numNodes);
      std::fill(tD.aPast.begin(),tD.aPast.end(),0);
      std::fill(tD.pPast.begin(),tD.pPast.end(),0);
    }
    else{
      tD.aPast = tD.a;
      tD.pPast = tD.p;
    }
    tD.a.clear();
    tD.p.clear();

    // reserve space
    sD.infected.reserve(numNodes);
    sD.notInfec.reserve(numNodes);
    sD.newInfec.reserve(numNodes);
    sD.status.reserve(numNodes);

    tD.a.reserve(numNodes);
    tD.p.reserve(numNodes);
    
    // begin filling in data
    for(j = 0; j < numNodes; ++j){
      hij = hist[i][j];
      
      sD.status.push_back(hij);
      if(hij < 2){ // not infected
	++sD.numNotInfec;
	sD.notInfec.push_back(j);
	
	tD.a.push_back(0);
	if(hij == 0){ // not treated
	  tD.p.push_back(0);
	}
	else{ // treated
	  tD.p.push_back(1);
	}
      }
      else{ // infected
	++sD.numInfected;
	sD.infected.push_back(j);

	if(++sD.timeInf[j] == 1){
	  sD.newInfec.push_back(j);
	}

	tD.p.push_back(0);
	if(hij == 2){ // not treated
	  tD.a.push_back(0);
	}
	else{ // treated
	  tD.a.push_back(1);
	}
      }
    }

    res.push_back(DataBundle(sD,tD,dD));
  }
  return res;
}
