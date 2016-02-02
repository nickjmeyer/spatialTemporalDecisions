#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef Model2GravityGDist M;
  typedef System<M,M> S;


  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);
  s.revert();

  Starts starts(1,s.fD.numNodes);
  s.reset(starts[0]);


  int i;
  for(i = 0; i < 15; ++i){
    s.nextPoint();
  }

  s.modelEst.estimateMle(s.sD,s.tD,s.fD,s.dD);

  std::cout << "gen: " << njm::toString(s.modelGen.getPar()) << std::endl;
  std::cout << "est: " << njm::toString(s.modelEst.getPar()) << std::endl;

  s.modelEst.setFisher(s.sD,s.tD,s.fD,s.dD);

  njm::sett.clean();
  return 0;
}
