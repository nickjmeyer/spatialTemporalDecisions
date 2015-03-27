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
  s.reset();

  omp_set_num_threads(4);

  std::cout << vr.run(s,nt,numReps,s.fD.finalT) << std::endl;
  
  // typedef GravityTimeInfExpCavesModel GM;
  // typedef GravityTimeInfExpCavesParam GP;
  // typedef GM EM;
  // typedef GP EP;

  // typedef System<GM,GP,EM,EP> S;

  // typedef ToyFeatures2<EM,EP> F;
  // typedef RankAgent<F,EM,EP> RA;

  // typedef FeaturesInt<F,EM,EP> FI;
  // typedef M2QOptim<S,RA,FI,EM,EP> OQ;

  // S s;

  // RA ra;
  // OQ oq;

  // s.modelGen.fitType = MLE;
  // s.modelEst.fitType = MLE;

  // s.paramEst_r = s.paramGen_r;
  // s.reset();

  // int i;
  // for(i = 0; i < 1; i++)
  //   njm::runif01();
  
  // int t;
  // for(t = 0; t < 10; ++t){
  //   if(t >= s.fD.trtStart){
  //     ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.paramEst);
  //   }
  //   s.nextPoint();

  // }

  // std::cout << "value: " << s.value() << std::endl;
  
  // oq.qEval.preCompData(s.sD,s.fD);

  // oq.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst);

  // oq.qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,ra);

  // oq.qEval.buildRD();

  // for(i = 0; i < 10; ++i)
  //   oq.qEval.solve();

  njm::sett.clean();
  return 0;
}
