#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesModel ME;

  typedef System<MG,ME> S;
  
  typedef ToyFeatures2<ME> F;

  // typedef NoTrt<ME> NT;
  typedef OsspAgent<ME> OA;

  typedef M1OsspOptim<S,OA,F,ME> OO;

  // typedef VanillaRunnerNS<S,NT> VNT;
  // typedef VanillaRunnerNS<S,RA> VRA;
  typedef OptimRunnerNS<S,OA,OO> OR;

  S s;

  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 2;
  Starts starts(numReps,s.fD.numNodes);
 
  // NT nt;
  OA oa;

  OO oo;

  // VNT vnt;
  // VRA vra;
  OR oor;


  // njm::timer.start("no trt");
  // vnt.run(s,nt,numReps,s.fD.finalT,starts);
  // njm::timer.stop("no trt");

  omp_set_num_threads(1);

  njm::timer.start("rank");
  oor.run(s,oa,oo,numReps,s.fD.finalT,starts);
  njm::timer.stop("rank");


  /////////////////////////////////////////////////////////

  // int i;
  // s.reset(starts[0]);
  // s.revert();
  // for(i = 0; i < 7; ++i){
  //   s.nextPoint();
  // }

  // std::cout << "prop inf: " << s.value() << std::endl;

  // OsspAgent<ME> oa;
  // M1OsspOptim<S,OsspAgent<ME>,F,ME>  oo;

  // oo.optim(s,oa);
  
  
  njm::sett.clean();
  return 0;
}
