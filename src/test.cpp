#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);



  typedef GravityModel MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  typedef ToyFeatures2<MG> F;
  typedef RankAgent<F,ME> RA;
  typedef PlainRunner<S,RA> PR;

  S s;
  RA ra;
  PR pr;

  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  s.modelEst_r = s.modelGen_r;
  

  int numReps = 1;
  Starts starts("startingLocations.txt");

  s.reset(starts[0]);
  
  njm::message(pr.run(s,ra,numReps,s.fD.finalT));
  // std::cout << "\n\nFUCKERS\n\n";
  // njm::message(pr.run(s,ra,1,s.fD.finalT));
  

  njm::sett.clean();
  return 0;
}
