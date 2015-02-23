#include "test3.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);


  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef System<GM,GP,GM,GP> S;

  S s("dataObs.txt");
  
  GravityMcmc mcmc;

  mcmc.load(s.sD.history,s.sD.status,s.fD);
  mcmc.sample(5000,1000);

  // {
  //   typedef GravityModel GM;
  //   typedef GravityParam GP;
  //   typedef System<GM,GP,GM,GP> S;

  //   S s;
  //   s.paramEst_r = s.paramGen_r;
  //   s.reset();

  //   for(int i = 0; i < 15; i++){
  //     s.updateStatus();
  //     s.nextPoint();
  //   }

  //   njm::message("Model: GRAVITY");
  //   njm::message("  Inf: " + njm::toString(s.value()));
  //   njm::message("Truth: " + njm::toString(s.paramGen.getPar()," ","\n"));
    
  //   s.modelEst.fitType = MLE;
  //   s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
  //   njm::message("  MLE: " + njm::toString(s.paramEst.getPar()," ","\n"));

  //   s.modelEst.fitType = MCMC;
  //   s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
  //   njm::message(" MCMc: " + njm::toString(s.paramEst.getPar()," ","\n"));
  // }


  // {
  //   typedef RangeModel GM;
  //   typedef RangeParam GP;
  //   typedef System<GM,GP,GM,GP> S;

  //   S s;
  //   s.paramGen_r.putPar({-3.0,300.0,1.0,4.0,4.0});
  //   s.paramEst_r = s.paramGen_r;
  //   s.reset();

  //   for(int i = 0; i < 15; i++){
  //     s.updateStatus();
  //     s.nextPoint();
  //   }

  //   njm::message("Model: RANGE");
  //   njm::message("  Inf: " + njm::toString(s.value()));
  //   njm::message("Truth: " + njm::toString(s.paramGen.getPar()," ","\n"));
    
  //   s.modelEst.fitType = MLE;
  //   s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
  //   njm::message("  MLE: " + njm::toString(s.paramEst.getPar()," ","\n"));

  //   s.modelEst.fitType = MCMC;
  //   s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
  //   njm::message(" MCMc: " + njm::toString(s.paramEst.getPar()," ","\n"));
  // }

  // {
  //   typedef CaveModel GM;
  //   typedef CaveParam GP;
  //   typedef System<GM,GP,GM,GP> S;

  //   S s;
  //   s.paramGen_r.putPar({-3.0,0.5,4.0,4.0});
  //   s.paramEst_r = s.paramGen_r;
  //   s.reset();

  //   for(int i = 0; i < 15; i++){
  //     s.updateStatus();
  //     s.nextPoint();
  //   }

  //   njm::message("Model: CAVE");
  //   njm::message("  Inf: " + njm::toString(s.value()));
  //   njm::message("Truth: " + njm::toString(s.paramGen.getPar()," ","\n"));
    
  //   s.modelEst.fitType = MLE;
  //   s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
  //   njm::message("  MLE: " + njm::toString(s.paramEst.getPar()," ","\n"));

  //   s.modelEst.fitType = MCMC;
  //   s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,s.paramEst);
  //   njm::message(" MCMc: " + njm::toString(s.paramEst.getPar()," ","\n"));
  // }
  
  
  
  njm::sett.clean();
  return 0;
}
