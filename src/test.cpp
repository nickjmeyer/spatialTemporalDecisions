#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  MultiModel mm;

  typedef GravityModel M0;
    
  mm.m.push_back(M0());

  Starts starts("startingLocations.txt");
  System<M0,M0> s;
  s.reset(starts[0]);
  s.revert();

  for(int i = 0; i < 10; ++i)
    s.nextPoint();

  mm.m.at(0).load(s.sD,s.tD,s.fD,s.dD);
  
  njm::sett.clean();
  return 0;
}
