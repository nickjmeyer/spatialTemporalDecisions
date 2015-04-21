#include "runM1MleMiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  
  typedef RadiusModel ME;

  typedef System<MG,ME> S;

  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures2<ME> F;
  // typedef RankAgent<F,ME> RA;
  typedef OsspAgent<ME> OA;
  
  // typedef M1SpOptim<S,RA,ME> SPO;
  typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  // typedef OptimRunner<S,RA,SPO> R_RA;
  typedef OptimRunner<S,OA,OSSPO> R_OA;


  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  MA ma;
  // RA ra;
  OA oa;

  // ra.tp.weights_r.zeros(ra.f.numFeatures);
  // ra.tp.weights_r(2) = 1;
  
  // SPO spo;
  // // no tuning for right now....
  // spo.tp.tune = 0;
  OSSPO osspo;

  R_MA r_ma;
  // R_RA r_ra;
  R_OA r_oa;
  

  njm::message("        Myopic: "
	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT,starts),
			       ""));
  njm::message("Priority Score: "
	       + njm::toString(r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts),
			       ""));

  return 0;
}

