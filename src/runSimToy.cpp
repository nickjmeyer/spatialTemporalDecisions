#include "runSimToy.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel MG;
  typedef GravityParam PG;
  
  // typedef GravityModel ME;
  // typedef GravityParam PE;
  typedef RangeModel ME;
  typedef RangeParam PE;

  typedef System<MG,PG,ME,PE> S;
  
  typedef ToyFeatures2<ME,PE> F;
  typedef FeaturesInt<F,ME,PE> FI;
  
  typedef NoTrt<ME,PE> AN;
  typedef ProximalAgent<ME,PE> AP;
  typedef MyopicAgent<ME,PE>  AM;
  typedef RankToyAgent<F,ME,PE> AR;

  typedef M1SgdOptim<S,AR,ME,PE> OM1_Sgd;
  typedef M2SaOptim<S,AR,FI,ME,PE> OM2_Sa;

  typedef VanillaRunner<S,AN> R_AN;
  typedef VanillaRunner<S,AP> R_AP;
  typedef FitOnlyRunner<S,AM> R_AM;
  typedef OptimRunner<S,AR,OM1_Sgd> R_AR_M1Sgd;
  typedef OptimRunner<S,AR,OM2_Sa> R_AR_M2Sa;

  // system
  S s;

  // agents
  AN an;
  AP ap;
  AM am;
  AR ar;

  // optim
  OM1_Sgd om1_sgd;
  OM2_Sa om2_sa;

  // runners
  R_AN r_an;
  R_AP r_ap;
  R_AM r_am;
  R_AR_M1Sgd r_ar_m1sgd;
  R_AR_M2Sa r_ar_m2sa;


  int mcReps=300,numPoints = s.fD.finalT;
  njm::message("No Trt");
  njm::message(r_an.run(s,an,mcReps,numPoints));
  
  njm::message("Proximal");
  njm::message(r_ap.run(s,ap,mcReps,numPoints));
  
  njm::message("Myopic");
  njm::message(r_am.run(s,am,mcReps,numPoints));
  
  njm::message("Rank M1");
  njm::message(r_ar_m1sgd.run(s,ar,om1_sgd,mcReps,numPoints));

  njm::message("Rank M2");
  njm::message(r_ar_m2sa.run(s,ar,om2_sa,mcReps,numPoints));


  // njm::sett.clean();
  
  return 0;
}
