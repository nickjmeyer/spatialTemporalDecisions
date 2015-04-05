#include "runM1MleWnsMiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  
  typedef RadiusModel ME;

  typedef System<MG,ME> S;

  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures2<ME> F;
  typedef RankAgent<F,ME> RA;

  typedef M1SpOptim<S,RA,ME> SPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts("startingLocations.txt");

  MA ma;
  RA ra;

  SPO spo;
  // no tuning for right now....
  spo.tp.tune = 0;

  R_MA r_ma;
  R_RA r_ra;
  

  njm::message("        Myopic: "
	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT,starts),
			       ""));
  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts),
			       ""));

  return 0;
}

