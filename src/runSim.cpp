#include "runSim.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  
  NoTrt<GravityModel,GravityParam> nA;
  ProximalAgent<GravityModel,GravityParam> pA;
  MyopicAgent<GravityModel,GravityParam> mA;
  RankAgent<GravityModel,GravityParam> rA;

  VanillaRunner<System,NoTrt,GravityModel,GravityParam> nR;
  VanillaRunner<System,ProximalAgent,GravityModel,GravityParam> pR;
  FitOnlyRunner<System,MyopicAgent,GravityModel,GravityParam> mR;
  FitOnlyRunner<System,RankAgent,GravityModel,GravityParam> rR;
  OptimRunner<System,RankToyAgent,GravityModel,GravityParam,M1SimpleOptim> rR1;
  OptimRunner<System,RankToyAgent,GravityModel,GravityParam,M1SgdOptim> rR2;

  M1SimpleOptim<System,RankToyAgent,GravityModel,GravityParam> m1Simple;
  M1SgdOptim<System,RankToyAgent,GravityModel,GravityParam> m1Sgd;
  
  int mcReps=120,numPoints = s.fD.finalT;
  njm::message("No Treatment");
  njm::message(nR.run(s,nA,mcReps,numPoints));
  njm::message("Proximal");
  njm::message(pR.run(s,pA,mcReps,numPoints));
  njm::message("Myopic");
  njm::message(mR.run(s,mA,mcReps,numPoints));
  njm::message("Rank");
  rA.tp.weights.zeros();
  rA.tp.weights(0)=1;
  njm::message(rR.run(s,rA,mcReps,numPoints));
  rA.tp.weights.ones();
  njm::message("Rank Simple");
  njm::message(rR1.run(s,rA,m1Simple,mcReps,numPoints));
  njm::message("Rank M1");
  njm::message(rR2.run(s,rA,m1Sgd,mcReps,numPoints));

  // njm::sett.clean();
  return 0;
}
