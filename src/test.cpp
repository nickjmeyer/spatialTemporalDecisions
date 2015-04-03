#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel M;
  typedef System<M,M> S;
  typedef NoTrt<M> N;
  typedef VanillaRunner<S,N> R;

  S s;
  N n;
  R r;

  int numReps = 500;
  Starts starts(numReps,s.fD.numNodes);
  
  std::cout << r.run(s,n,500,s.fD.finalT,starts) << std::endl;

  
  njm::sett.clean();
  return 0;
}
