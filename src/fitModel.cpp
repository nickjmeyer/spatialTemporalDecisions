#include "fitModel.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    njm::message("Gravity Model without time infected");
    
    typedef GravityModel GM;
    typedef GravityParam GP;
    typedef GravityModel EM;
    typedef GravityParam EP;

    typedef System<GM,GP,EM,EP> S;
  
    S s("dataObs.txt");
  
    s.paramEst_r = s.paramGen_r;
    s.reset();

    s.modelGen.fitType = MLE;
    s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramGen);

    njm::message(" MLE: " + njm::toString(s.paramGen.getPar()," "));

    s.modelGen.fitType = MCMC;
    s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramGen);

    njm::message("MCMC: " + njm::toString(s.paramGen.getPar()," "));

    s.paramGen.save();
  }


  {
    njm::message("Gravity Model with time infected");
    
    typedef GravityTimeInfModel GM;
    typedef GravityTimeInfParam GP;
    typedef GravityTimeInfModel EM;
    typedef GravityTimeInfParam EP;

    typedef System<GM,GP,EM,EP> S;
  
    S s("dataObs.txt");
  
    s.paramEst_r = s.paramGen_r;
    s.reset();

    s.modelGen.fitType = MLE;
    s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramGen);

    njm::message(" MLE: " + njm::toString(s.paramGen.getPar()," "));

    s.modelGen.fitType = MCMC;
    s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramGen);

    njm::message("MCMC: " + njm::toString(s.paramGen.getPar()," "));

    s.paramGen.save();
  }


  njm::sett.clean();
  return 0;
}
