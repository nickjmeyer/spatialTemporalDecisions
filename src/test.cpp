#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCaves MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  typedef ProximalAgent<ME> PA;

  S s;
  PA pa;

  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  s.modelEst_r = s.modelGen_r;
  s.revert();

  Starts starts("startingLocations.txt");
  s.reset(starts[0]);
  s.revert();

  int i;
  for(i = 0; i < 2; ++i){
    if(i >= s.fD.trtStart)
      pa.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
    s.updateStatus();

    s.nextPoint();
  }

  s.modelGen.setType(MLES);
  s.modelGen.setFisher(s.sD,s.tD,s.fD,s.dD);
  for(i = 0; i < 20; ++i){
    s.modelGen.sample();
    std::cout << njm::toString(s.modelGen.getPar()," ","\n");
  }


  njm::sett.clean();
  return 0;
}
