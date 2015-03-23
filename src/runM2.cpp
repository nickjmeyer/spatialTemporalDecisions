#include "runM2.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef MG ME;
  typedef PG PE;

  typedef System<MG,PG,ME,PE> S;

  typedef ToyFeatures2<ME,PE> F;
  typedef FeaturesInt<F,ME,PE> FI;
  typedef RankAgent<F,ME,PE> RA;

  typedef M2QOptim<S,RA,FI,ME,PE> SPO;

  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;

  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  RA ra;

  SPO spo;

  R_RA r_ra;
  

  int numReps = 100;

  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT),""));

  return 0;
}

