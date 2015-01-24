#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam_r = s.genParam_r;

  ProximalAgent<GravityModel,GravityParam> pA;
  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
	       GravityModel,GravityParam> rA;
  
  PlainRunner<System<GravityModel,GravityParam>,
  	      ProximalAgent<GravityModel,GravityParam> > nR;

  int i;
  for(i=0; i<12; i++){
    if(i>=s.fD.trtStart)
      pA.applyTrt(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);

    s.updateStatus();
    s.nextPoint();
  }

  M2NmOptim<System<GravityModel,GravityParam>,
	    RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
			 GravityModel,GravityParam>,
	    ToyFeatures2<GravityModel,GravityParam>,
	    GravityModel,GravityParam> m2;

  m2.qEval.preCompData(s.sD,s.fD);
  m2.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  m2.qEval.solve();

  njm::message(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA));

  m2.optim(s,rA);

  m2.qEval.preCompData(s.sD,s.fD);
  m2.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  m2.qEval.solve();

  njm::message(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA));


  njm::sett.clean();
  return 0;
}
