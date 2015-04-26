#include "runM1MleMiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCaves MG;
  
  typedef ModelRadius ME;

  typedef System<MG,ME> S;

  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures3<ME> F;
  typedef OsspAgent<ME> OA;
  
  typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,OA,OSSPO> R_OA;


  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  MA ma;
  OA oa;

  OSSPO osspo;

  R_MA r_ma;
  R_OA r_oa;
  

  njm::message("        Myopic: "
	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT,starts),
			       ""));
  njm::message("Priority Score: "
	       + njm::toString(r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts),
			       ""));

  return 0;
}

