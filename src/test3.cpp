#include "test3.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef GravityModel GM;
    typedef GravityParam GP;
    typedef System<GM,GP,GM,GP> S;

    S s;
    s.paramEst_r = s.paramGen_r;
    s.reset();

    for(int i = 0; i < 15; i++){
      s.updateStatus();
      s.nextPoint();
    }

    njm::message("Truth: " + njm::toString(s.paramGen.getPar()," ","\n"));
    
    s.modelEst.fitType = MLE;
    s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
    njm::message("  MLE: " + njm::toString(s.paramEst.getPar()," ","\n"));

    s.modelEst.fitType = MCMC;
    s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
    njm::message(" MCMc: " + njm::toString(s.paramEst.getPar()," ","\n"));
  }
  
  njm::sett.clean();
  return 0;
}
