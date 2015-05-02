#include "runM1MlesWnsMiss.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCaves MG;
  
  typedef ModelDist ME;

  typedef System<MG,ME> S;

  typedef MyopicAgent<ME> MA;
  
  typedef WnsFeatures1<ME> F;
  typedef RankAgent<F,ME> RA;
  // typedef OsspAgent<ME> OA;

  typedef M1SpOptim<S,RA,ME> SPO;
  // typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;
  // typedef OptimRunner<S,OA,OSSPO> R_OA;


  // S s;
  S s("obsData.txt");
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
  Starts starts("startingLocations.txt");

  MA ma;
  RA ra;
  // OA oa;

  SPO spo;
  // OSSPO osspo;

  R_MA r_ma;
  R_RA r_ra;
  // R_OA r_oa;


  RunStats rs;

  std::vector<double> times = {15,25};
  std::vector<double>::const_iterator it,beg,end;
  beg = times.begin();
  end = times.end();
    
  for(it = beg; it != end; ++it){
    s.fD.finalT = s.sD.time + *it;

    ma.name = "myopic_obs" + njm::toString(*it,"",0,0);
    rs = r_ma.run(s,ma,numReps,s.fD.finalT,starts);
    njm::message("         Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    ra.name = "rank_obs" + njm::toString(*it,"",0,0);
    rs = r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts);
    njm::message("  Policy Search: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");
    
    // rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
    // njm::message("One Step Polish: "
    // 	       + njm::toString(rs.smean(),"")
    // 	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  }
  return 0;
}

