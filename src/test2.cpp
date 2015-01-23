#include "test2.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s0,s1;
  s0.estParam_r = s0.genParam_r;
  s1.estParam_r = s1.genParam_r;
  s0.reset();
  s1.reset();

  PlainRunner<System,RankToyAgent,GravityModel,GravityParam> pR;

  RankToyAgent<GravityModel,GravityParam> rA;

  M1NmOptim <System,RankToyAgent,GravityModel,GravityParam> qO;


  int t,T=10;
  for(t=0; t<T; t++){
    
    if(t>=s0.fD.trtStart)
      rA.applyTrt(s0.sD,s0.tD,s0.fD,s0.dD,s0.model,s0.genParam);
    if(t>=s1.fD.trtStart)
      rA.applyTrt(s1.sD,s1.tD,s1.fD,s1.dD,s1.model,s1.genParam);
    
    s0.updateStatus();
    s1.updateStatus();
    
    s0.nextPoint();
    s1.nextPoint();
  }

  
  qO.qEval.preCompData(s0.sD,s0.fD);
  qO.qEval.bellResFixData(s0.sD,s0.tD,s0.fD,s0.dD,s0.model,s0.genParam);
  qO.qEval.bellResPolData(s0.sD.time,s0.fD,s0.model,s0.genParam,rA);
  qO.qEval.solve();
  std::cout << qO.qEval.qFn(s0.sD,s0.tD,s0.fD,s0.dD,s0.model,s0.genParam,rA)
  	    << " || " << qO.qEval.bellRes() << std::endl;

  qO.qEval.preCompData(s1.sD,s1.fD);
  qO.qEval.bellResFixData(s1.sD,s1.tD,s1.fD,s1.dD,s1.model,s1.genParam);
  qO.qEval.bellResPolData(s1.sD.time,s1.fD,s1.model,s1.genParam,rA);
  qO.qEval.solve();
  std::cout << qO.qEval.qFn(s1.sD,s1.tD,s1.fD,s1.dD,s1.model,s1.genParam,rA)
  	    << " || " << qO.qEval.bellRes() << std::endl;

  qO.qEval.preCompData(s0.sD,s0.fD);
  qO.qEval.bellResFixData(s0.sD,s0.tD,s0.fD,s0.dD,s0.model,s0.genParam);
  qO.qEval.bellResPolData(s0.sD.time,s0.fD,s0.model,s0.genParam,rA);
  qO.qEval.solve();
  std::cout << qO.qEval.qFn(s0.sD,s0.tD,s0.fD,s0.dD,s0.model,s0.genParam,rA)
  	    << " || " << qO.qEval.bellRes() << std::endl;

  
  njm::sett.clean();
  return 0;
}
