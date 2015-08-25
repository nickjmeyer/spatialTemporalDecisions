#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef Model2GravityGDist M;
  typedef System<M,M> S;


  S s;
  Starts starts(1,s.fD.numNodes);
  s.modelEst_r = s.modelGen_r;
  s.reset(starts[0]);


  int i;
  for(i = 0; i < 2; ++i){
    s.nextPoint();
  }

  std::cout << s.value() << std::endl;

  RankAgent<ToyFeatures5<M>,M> ra;

  // ra.tp.jitterScale = -1;
  // ra.tp.shuffle = false;

  int j;
  for(i = 0; i < 10; ++i){
    std::fill(s.tD.a.begin(),s.tD.a.end(),0);
    std::fill(s.tD.p.begin(),s.tD.p.end(),0);
    ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen);

    for(j = 0; j < s.fD.numNodes; ++j){
      if(s.tD.a.at(j) == 1 || s.tD.p.at(j) == 1){
	std::cout << j << " ";
      }
    }
    std::cout << std::endl;
  }


  njm::sett.clean();
  return 0;
}
