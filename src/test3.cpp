#include "test3.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef System<GM,GP,GM,GP> S;

  S s("dataObs.txt");

  GravityMcmc mcmc;

  mcmc.load(s.sD.history,s.sD.status,s.fD);
  mcmc.sample(100,50);
  

  njm::sett.clean();
  return 0;
}
