#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef GravityModel EM;
  typedef GravityParam EP;

  typedef System<GM,GP,EM,EP> S;

  typedef ToyFeatures2<EM,EP> F;

  typedef RankToyAgent<F,EM,EP> RA;

  typedef PlainRunner<S,RA> PR;
  
  S s;
  RA ra;
  PR pr;
  
  s.paramEst_r = s.paramGen_r;
  s.reset();

  pr.run(s,ra,10,15);


  njm::sett.clean();
  return 0;
}
