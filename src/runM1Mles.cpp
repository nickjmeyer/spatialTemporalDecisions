#include "runM1Mles.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDistTrendPowCon MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef NoTrt<ME> NT;
  typedef ProximalGDistAgent<ME> PA;
  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures5<ME> F5;
  typedef ToyFeatures6<ME> F6;
  typedef RankAgent<F5,ME> RA5;
  typedef RankAgent<F6,ME> RA6;
  typedef OsspAgent<ME> OA;

  typedef M1SpOptim<S,RA5,ME> SPO5;
  typedef M1OsspOptim<S,OA,F5,ME> OSSPO5;
  typedef M1SpOptim<S,RA6,ME> SPO6;
  typedef M1OsspOptim<S,OA,F6,ME> OSSPO6;
  
  typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunner<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA5,SPO5> R_RA5;
  typedef OptimRunner<S,OA,OSSPO5> R_OA5;
  typedef OptimRunner<S,RA6,SPO6> R_RA6;
  typedef OptimRunner<S,OA,OSSPO6> R_OA6;


  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  NT nt;
  PA pa;
  MA ma;
  RA5 ra5;
  RA6 ra6;
  OA oa;

  ra5.name = "rank_5";
  ra6.name = "rank_6";
  

  SPO5 spo5;
  OSSPO5 osspo5;
  SPO6 spo6;
  OSSPO6 osspo6;

  osspo5.name = "M1Ossp_5";
  osspo6.name = "M1Ossp_6";

  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_RA5 r_ra5;
  R_RA6 r_ra6;
  R_OA5 r_oa5;
  R_OA6 r_oa6;


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

  rs = r_ra6.run(s,ra6,spo6,numReps,s.fD.finalT,starts);
  njm::message("  Policy Search: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  std::vector<int> N = {100, 10000};
  std::vector<double> JS = {4.0, 1.0};
  std::vector<int>::const_iterator itN;
  std::vector<double>::const_iterator itJS;

  for(itN = N.begin(); itN != N.end(); ++itN){
    for(itJS = JS.begin(); itJS != JS.end(); ++itJS){
  
      osspo5.tp.N = *itN;
      osspo5.tp.jitterScale = *itJS;
      osspo5.name = "M1Ossp_" + njm::toString(osspo5.tp.N,"",0,0)
  	+ "_" + njm::toString(osspo5.tp.jitterScale,"",0,0);
      rs = r_oa5.run(s,oa,osspo5,numReps,s.fD.finalT,starts);
      njm::message("One Step Polish: "
  		   + njm::toString(rs.smean(),"")
  		   + "  (" + njm::toString(rs.seMean(),"") + ")");
  
      osspo6.tp.N = *itN;
      osspo6.tp.jitterScale = *itJS;
      osspo6.name = "M1Ossp_" + njm::toString(osspo6.tp.N,"",0,0)
  	+ "_" + njm::toString(osspo6.tp.jitterScale,"",0,0);
      rs = r_oa6.run(s,oa,osspo6,numReps,s.fD.finalT,starts);
      njm::message("One Step Polish: "
  		   + njm::toString(rs.smean(),"")
  		   + "  (" + njm::toString(rs.seMean(),"") + ")");
    }
  }
  
  return 0;
}

