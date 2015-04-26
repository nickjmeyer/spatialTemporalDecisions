#include "runM1Mle.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCaves MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef NoTrt<ME> NT;
  typedef ProximalAgent<ME> PA;
  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures3<ME> F;
  typedef OsspAgent<ME> OA;

  typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunner<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,OA,OSSPO> R_OA;


  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  NT nt;
  PA pa;
  MA ma;
  OA oa;

  OSSPO osspo;

  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_OA r_oa;
  

    
  // njm::message("  No treatment: "
  // 	       + njm::toString(r_nt.run(s,nt,numReps,s.fD.finalT,starts),
  // 			       ""));
  // njm::message("      Proximal: "
  // 	       + njm::toString(r_pa.run(s,pa,numReps,s.fD.finalT,starts),
  // 			       ""));
  // njm::message("        Myopic: "
  // 	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT,starts),
  // 			       ""));
  njm::message("Priority Score: "
  	       + njm::toString(r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts),
  			       ""));
  return 0;
}

