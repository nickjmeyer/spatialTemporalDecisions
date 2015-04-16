#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);



  typedef GravityTimeInfExpCavesModel MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  typedef ToyFeatures2<MG> F;
  typedef RankAgent<F,ME> RA;
  typedef PlainRunner<S,RA> PR;

  S s;
  RA ra;
  PR pr;


  int numReps = 1000;
  Starts starts("startingLocations.txt");

  s.reset(starts[0]);

  pr.run(s,ra,numReps,s.fD.finalT);


  njm::sett.clean();
  return 0;
}
