#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<EbolaModel,EbolaParam> s;
  s.estParam_r = s.genParam_r;
  NoTrt<EbolaModel,EbolaParam> nA;
  
  PlainRunner<System<EbolaModel,EbolaParam>,
  	      NoTrt<EbolaModel,EbolaParam> > nR;

  njm::message(s.fD.trtStart);
  njm::message(s.fD.period);
  
  int N=300;
  njm::message(nR.run(s,nA,N,250));

  njm::sett.clean();
  return 0;
}
