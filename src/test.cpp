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

  typedef M1SgdOptim<S,RA,EM,EP> M1;

  typedef OptimRunner<S,RA,M1> OR;
  
  S s;
  RA ra;
  PR pr;
  M1 m1;
  OR o1;
  
  s.paramEst_r = s.paramGen_r;
  s.reset();

  njm::message(pr.run(s,ra,600,15));
  njm::message(o1.run(s,ra,m1,100,15));


  njm::sett.clean();
  return 0;
}
