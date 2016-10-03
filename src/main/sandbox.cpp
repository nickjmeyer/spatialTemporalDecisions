#include <gflags/gflags.h>
#include <glog/logging.h>
#include "tuneGenWNS.hpp"

DEFINE_string(srcDir,"","Path to source directory");

int main(int argc, char ** argv){
  ::google::InitGoogleLogging(argv[0]);
  ::gflags::ParseCommandLineFlags(&argc,&argv,true);
  njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

  const bool edgeToEdge = false;

  typedef Model2GravityEDist MG;

  typedef MG ME;

  typedef System<MG,ME> S;

  typedef WnsFeatures3<ME> F;
  typedef RankAgent<F,ME> RA;

  typedef VanillaRunnerNS<S,RA> RR;

  S s("obsData.txt");
  s.setEdgeToEdge(edgeToEdge);
  s.modelEst_r = s.modelGen_r;
  s.revert();

  int numReps = 500;
  Starts starts("startingLocations.txt");


  RA ra;
  RR rr;
  ra.setEdgeToEdge(edgeToEdge);
  // ra.reset();

  double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).smean();
  std::cout << " valRA: " + njm::toString(valRA,"") << std::endl;

  njm::sett.clean();

  return 0;
}
