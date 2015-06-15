#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDistTrendPowCon MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef ToyFeatures5<ME> F5;

  // typedef ProximalGDistAgent<ME> PA;
  typedef RankAgent<F5,ME> RA5;

  // typedef NullOptim<S,PA,ME> NO;
  typedef M1SpOptim<S,RA5,ME> RSA;

  // typedef IncremAgent<ME,PA,NO> IA;
  typedef IncremAgent<ME,RA5,RSA> IA;

  typedef FitOnlyRunner<S,IA> R_IA;
  

  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
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
