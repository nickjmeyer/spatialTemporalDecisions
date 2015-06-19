#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDistTrendPowCon MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef NoTrt<ME> NT;
  typedef ProximalGDistAgent<ME> PA;
  typedef MyopicAgent<ME> MA;

  typedef ToyFeatures5<ME> F5;
  typedef RankAgent<F5,ME> RA5;

  typedef M1SpOptim<S,RA5,ME> SPO5;
  
  typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunner<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA5,SPO5> R_RA5;

  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  NT nt;
  PA pa;
  MA ma;
  RA5 ra5;

  ra5.name = "rank_5";
  ra5.tp.jitterScale = -1;

  SPO5 spo5;
  spo5.tp.fixSample = 1;
  
  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_RA5 r_ra5;

  RunStats rs;

  rs = r_nt.run(s,nt,numReps,s.fD.finalT,starts);
  njm::message("   No treatment: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_pa.run(s,pa,numReps,s.fD.finalT,starts);
  njm::message("       Proximal: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_ma.run(s,ma,numReps,s.fD.finalT,starts);
  njm::message("         Myopic: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_ra5.run(s,ra5,spo5,numReps,s.fD.finalT,starts);
  njm::message("  Policy Search: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  return 0;
}
