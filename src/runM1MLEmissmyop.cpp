#include "runM1MLEmissmyop.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef RadiusModel ME;
  typedef RadiusParam PE;

  typedef System<MG,PG,ME,PE> S;

  typedef MyopicAgent<ME,PE> MA;
  
  typedef FitOnlyRunner<S,MA> R_MA;


  S s;
  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  MA ma;

  R_MA r_ma;

  int numReps = 100;
  

  njm::message("        Myopic: "
	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT),""));

  return 0;
}

