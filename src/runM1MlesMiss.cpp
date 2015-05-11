#include "runM1MlesMiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);


  typedef ModelTimeExpCavesGDist MG;
  
  typedef ModelTimeExpCavesEDist ME;

  typedef System<MG,ME> S;

  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures4<ME> F;
  typedef RankAgent<F,ME> RA;
  typedef OsspAgent<ME> OA;

  typedef M1SpOptim<S,RA,ME> SPO;
  typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;
  typedef OptimRunner<S,OA,OSSPO> R_OA;


  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  MA ma;
  RA ra;
  OA oa;

  SPO spo;
  OSSPO osspo;

  R_MA r_ma;
  R_RA r_ra;
  R_OA r_oa;


  RunStats rs;

  rs = r_ma.run(s,ma,numReps,s.fD.finalT,starts);
  njm::message("         Myopic: "
  	       + njm::toString(rs.smean(),"")
	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts);
  njm::message("  Policy Search: "
  	       + njm::toString(rs.smean(),"")
	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
  njm::message("One Step Polish: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  return 0;
}

