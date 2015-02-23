#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef GravityModel GM;
    typedef GravityParam GP;
    typedef GravityModel EM;
    typedef GravityParam EP;

    typedef System<GM,GP,EM,EP> S;
  
    S s;
  
    s.paramEst_r = s.paramGen_r;
    s.reset();

    int i;
    for(i = 0; i < 15; i++){
      s.updateStatus();
      s.nextPoint();
    }

    s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);

    njm::message(s.paramGen.getPar());
  }


  njm::sett.clean();
  return 0;
}
