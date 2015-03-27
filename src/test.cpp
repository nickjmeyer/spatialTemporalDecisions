#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel GM;
  typedef GravityTimeInfExpCavesParam GP;
  typedef GM EM;
  typedef GP EP;

  typedef System<GM,GP,EM,EP> S;

  typedef NoTrt<EM,EP> NT;

  typedef VanillaRunner<S,NT> VR;

  S s;


  int numReps = 500;
  Starts starts(numReps,s.fD.numNodes);
  // Starts starts("startingLocations.txt");

  NT nt;
  VR vr;

  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  s.paramEst_r = s.paramGen_r;
  s.revert();

  omp_set_num_threads(1);

  std::cout << vr.run(s,nt,numReps,s.fD.finalT,starts) << std::endl;
  std::cout << vr.run(s,nt,numReps,s.fD.finalT,starts) << std::endl;
  std::cout << vr.run(s,nt,numReps,s.fD.finalT,starts) << std::endl;
  
  njm::sett.clean();
  return 0;
}
