#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam_r = s.genParam_r;
  s.reset();

  ProximalAgent<GravityModel,GravityParam> pA;
  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
	       GravityModel,GravityParam> rA;

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
  	    FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
  			GravityModel,GravityParam>,
  	    GravityModel,GravityParam> m2;

  // std::string D,P,R,B,Q;

  // m2.qEval.preCompData(s.sD,s.fD);
  // m2.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  // m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  // m2.qEval.solve();
  // D = njm::toString(m2.qEval.D.sum(),"",8,4);
  // P = njm::toString(m2.qEval.P.sum(),"",8,4);
  // R = njm::toString(m2.qEval.R.sum(),"",8,4);
  // B = njm::toString(m2.qEval.beta.sum(),"",8,4);
  // Q = njm::toString(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA),
  // 		    "",8,4);
  // njm::message("one");
  // njm::message("D: " + D + " -- P: " + P + " -- R: " + R + " -- B: " + B +
  // 	       " -- Q: " + Q);
  

  // for(i=0; i<rA.f.numFeatures; i++)
  //   rA.tp.weights.at(i)+= njm::rnorm01();

  // m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  // m2.qEval.solve();
  // D = njm::toString(m2.qEval.D.sum(),"",8,4);
  // P = njm::toString(m2.qEval.P.sum(),"",8,4);
  // R = njm::toString(m2.qEval.R.sum(),"",8,4);
  // B = njm::toString(m2.qEval.beta.sum(),"",8,4);
  // Q = njm::toString(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA),
  // 		    "",8,4);
  // njm::message("two");
  // njm::message("D: " + D + " -- P: " + P + " -- R: " + R + " -- B: " + B +
  // 	       " -- Q: " + Q);


  // m2.qEval.preCompData(s.sD,s.fD);
  // m2.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  // m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  // m2.qEval.solve();
  // D = njm::toString(m2.qEval.D.sum(),"",8,4);
  // P = njm::toString(m2.qEval.P.sum(),"",8,4);
  // R = njm::toString(m2.qEval.R.sum(),"",8,4);
  // B = njm::toString(m2.qEval.beta.sum(),"",8,4);
  // Q = njm::toString(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA),
  // 		    "",8,4);
  // njm::message("two");
  // njm::message("D: " + D + " -- P: " + P + " -- R: " + R + " -- B: " + B +
  // 	       " -- Q: " + Q);


  // m2.qEval.solve();
  // D = njm::toString(m2.qEval.D.sum(),"",8,4);
  // P = njm::toString(m2.qEval.P.sum(),"",8,4);
  // R = njm::toString(m2.qEval.R.sum(),"",8,4);
  // B = njm::toString(m2.qEval.beta.sum(),"",8,4);
  // Q = njm::toString(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA),
  // 		    "",8,4);
  // njm::message("two");
  // njm::message("D: " + D + " -- P: " + P + " -- R: " + R + " -- B: " + B +
  // 	       " -- Q: " + Q);
  
  // m2.qEval.solve();
  // D = njm::toString(m2.qEval.D.sum(),"",8,4);
  // P = njm::toString(m2.qEval.P.sum(),"",8,4);
  // R = njm::toString(m2.qEval.R.sum(),"",8,4);
  // B = njm::toString(m2.qEval.beta.sum(),"",8,4);
  // Q = njm::toString(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA),
  // 		    "",8,4);
  // njm::message("two");
  // njm::message("D: " + D + " -- P: " + P + " -- R: " + R + " -- B: " + B +
  // 	       " -- Q: " + Q);

  // m2.qEval.solve();
  // D = njm::toString(m2.qEval.D.sum(),"",8,4);
  // P = njm::toString(m2.qEval.P.sum(),"",8,4);
  // R = njm::toString(m2.qEval.R.sum(),"",8,4);
  // B = njm::toString(m2.qEval.beta.sum(),"",8,4);
  // Q = njm::toString(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA),
  // 		    "",8,4);
  // njm::message("two");
  // njm::message("D: " + D + " -- P: " + P + " -- R: " + R + " -- B: " + B +
  // 	       " -- Q: " + Q);

  //////////////////////////////////////////////////////////////////////////////

  m2.qEval.preCompData(s.sD,s.fD);
  m2.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  m2.qEval.solve();

  
  njm::message("before");
  njm::message(rA.tp.getPar());
  njm::message(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA));

  m2.optim(s,rA);

  m2.qEval.preCompData(s.sD,s.fD);
  m2.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  m2.qEval.bellResPolData(s.sD.time,s.fD,s.model,s.estParam,rA);
  m2.qEval.solve();

  njm::message("after");
  njm::message(rA.tp.getPar());
  njm::message(m2.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam,rA));


  njm::sett.clean();
  return 0;
}
