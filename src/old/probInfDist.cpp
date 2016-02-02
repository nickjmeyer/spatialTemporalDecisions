#include "probInfDist.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDist MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef ProximalGDistAgent<ME> PA;

  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 50;
  Starts starts(numReps,s.fD.numNodes);

  
  PA pa;

  int r,t,i;
  std::vector<double> probs;
  std::vector<int> trt;
  for(r = 0; r < numReps; ++r){
    njm::resetSeed(r);
    s.reset(starts[r]);
    for(t = s.sD.time; t < numPoints; ++t){
      if(t>=s.fD.trtStart){
	pa.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
      }

      s.updateStatus();

      trt.clear();
      probs.clear();
      for(i = 0; i < s.fD.numNodes; ++i){
	
      }

      s.nextPoint();
    }
  }


  njm::sett.clean();
  return 0;
}
