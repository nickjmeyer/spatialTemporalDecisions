#include "runM1MlesMiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);


  typedef ModelTimeExpCavesGPowGDistTrendPowCon MG;

  typedef ModelGDist ME;

  typedef System<MG,ME> S;

  typedef MyopicAgent<ME> MA;

  typedef ToyFeatures5<ME> F;
  typedef RankAgent<F,ME> RA;

  typedef M1SpOptim<S,RA,ME> SPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  MA ma;
  RA ra;

  ra.tp.jitterScale = -1;

  SPO spo;
  spo.tp.fixSample = 1;

  R_MA r_ma;
  R_RA r_ra;


  RunStats rs;

  rs = r_ma.run(s,ma,numReps,s.fD.finalT,starts);
  njm::message("         Myopic: "
  	       + njm::toString(rs.smean(),"")
	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  rs = r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts);
  njm::message("  Policy Search: "
  	       + njm::toString(rs.smean(),"")
	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  return 0;
}
