#include "runSimToy.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel MG;
  typedef GravityParam PG;
  typedef GravityModel ME;
  typedef GravityParam PE;

  typedef System<MG,PG,ME,PE> S;
  
  typedef ToyFeatures2<ME,PE> F;
  
  typedef NoTrt<ME,PE> AN;
  typedef ProximalAgent<ME,PE> AP;
  typedef MyopicAgent<ME,PE>  AM;
  typedef RankToyAgent<F,ME,PE> AR;

  typedef M1SimpleOptim<S,AR,ME,PE> OM1_Simple;
  typedef M1SgdOptim<S,AR,ME,PE> OM1_Sgd;
  typedef M1HybridOptim<S,AR,ME,PE> OM1_Hybrid;

  typedef VanillaRunner<S,AN> R_AN;
  typedef VanillaRunner<S,AP> R_AP;
  typedef FitOnlyRunner<S,AM> R_AM;
  typedef FitOnlyRunner<S,AR> R_AR_Fix;
  typedef OptimRunner<S,AR,OM1_Simple> R_AR_Simple;
  typedef OptimRunner<S,AR,OM1_Sgd> R_AR_Sgd;
  typedef OptimRunner<S,AR,OM1_Hybrid> R_AR_Hybrid;

  // system
  S s;

  // agents
  AN an;
  AP ap;
  AM am;
  AR ar;

  // optim
  OM1_Simple om1_simple;
  OM1_Sgd om1_sgd;
  OM1_Hybrid om1_hybrid;

  // runners
  R_AN r_an;
  R_AP r_ap;
  R_AM r_am;
  R_AR_Fix r_ar_fix;
  R_AR_Simple r_ar_simple;
  R_AR_Sgd r_ar_sgd;
  R_AR_Hybrid r_ar_hybrid;
  

    
  

  int mcReps=300,numPoints = s.fD.finalT;
  njm::message("No Treatment");
  njm::message(r_an.run(s,an,mcReps,numPoints));
  
  njm::message("Proximal");
  njm::message(r_ap.run(s,ap,mcReps,numPoints));
  
  njm::message("Myopic");
  njm::message(r_am.run(s,am,mcReps,numPoints));
  
  njm::message("Rank");
  njm::message(r_ar_fix.run(s,ar,mcReps,numPoints));
  
  njm::message("Rank Simple");
  njm::message(r_ar_simple.run(s,ar,om1_simple,mcReps,numPoints));

  njm::message("Rank Sgd");
  njm::message(r_ar_sgd.run(s,ar,om1_sgd,mcReps,numPoints));
  
  njm::message("Rank Hybrid");
  njm::message(r_ar_hybrid.run(s,ar,om1_hybrid,mcReps,numPoints));


  // njm::sett.clean();
  
  return 0;
}
