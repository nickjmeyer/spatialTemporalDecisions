#include "test3.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef System<GM,GP,GM,GP> System_;
  typedef ToyFeatures2<GM,GP> F;
  typedef FeaturesInt<F,GM,GP> FI;
  typedef RankToyAgent<F,GM,GP> Agent_;
  // typedef NoTrt<GM,GP> Agent_;
  typedef AnchorMan<System_,Agent_,FI,GM,GP> Optim_;
  typedef OptimRunner<System_,Agent_,Optim_> Runner_;
  // typedef PlainRunner<System_,Agent_> Runner_;

  System_ s;
  Agent_ a;
  Runner_ r;
  Optim_ o;

  njm::message(r.run(s,a,o,16,15));
  // njm::message(r.run(s,a,1000,15))
  
  njm::sett.clean();
  return 0;
}
