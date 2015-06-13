#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDistTrendPowCon MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef ProximalGDistAgent<ME> PA;

  typedef NullOptim<S,PA,ME> NO;

  typedef IncremAgent<ME,PA,NO> IA;

  typedef FitOnlyRunner<S,IA> R_IA;
  

  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 1;
  Starts starts(numReps,s.fD.numNodes);

  IA ia;
  
  R_IA r_ia;

  RunStats rs;

  rs = r_ia.run(s,ia,numReps,s.fD.finalT,starts);
  njm::message("    IncremAgent: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  njm::sett.clean();
  return 0;
}
