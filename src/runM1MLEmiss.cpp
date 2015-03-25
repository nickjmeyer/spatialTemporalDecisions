#include "runM1MLEmiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef RadiusModel ME;
  typedef RadiusParam PE;

  typedef System<MG,PG,ME,PE> S;

  typedef MyopicAgent<ME,PE> MA;
  
  typedef ToyFeatures2<ME,PE> F;
  typedef RankAgent<F,ME,PE> RA;

  typedef M1SpOptim<S,RA,ME,PE> SPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  MA ma;
  RA ra;

  ra.tp.weights_r.zeros(ra.f.numFeatures);
  ra.tp.weights_r(2) = 1;
  
  SPO spo;
  // no tuning for right now....
  spo.tp.tune = 0;

  R_MA r_ma;
  R_RA r_ra;
  

  int numReps = 100;
  

  njm::message("        Myopic: "
	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT),""));
  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT),""));

  return 0;
}

