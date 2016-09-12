#include <gflags/gflags.h>
#include <glog/logging.h>
#include "runM1MlesWns.hpp"

DEFINE_string(srcDir,"","Path to source directory");
// DEFINE_bool(edgeToEdge,false,"Edge to edge transmission");

int main(int argc, char ** argv){
  ::google::InitGoogleLogging(argv[0]);
  ::google::ParseCommandLineFlags(&argc,&argv,true);
  njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

  // typedef ModelTimeExpCavesGDistTrendPowCon MG;

  typedef Model2GravityGDist MG;

  typedef MG ME;

  typedef System<MG,ME> S;

  typedef NoTrt<ME> NT;
  typedef ProximalAgent<ME> PA;
  typedef MyopicAgent<ME> MA;

  typedef WnsFeatures3<ME> F;
  typedef RankAgent<F,ME> RA;

  typedef M1SpOptim<S,RA,ME> SPO;

  typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunner<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;

  const bool edgeToEdge = false;


  // S s;
  S s("obsData.txt");
  s.setEdgeToEdge(edgeToEdge);
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 100;
  Starts starts("startingLocations.txt");

  NT nt;
  PA pa;
  MA ma;
  RA ra;

  pa.setEdgeToEdge(edgeToEdge);

  ra.tp.jitterScale = -1;
  ra.setEdgeToEdge(edgeToEdge);

  SPO spo;
  spo.tp.fixSample = 1;

  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_RA r_ra;


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

  rs = r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts);
  njm::message("  Policy Search: "
	       + njm::toString(rs.smean(),"")
	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  return 0;
}
