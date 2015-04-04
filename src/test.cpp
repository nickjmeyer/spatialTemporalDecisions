#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  
  typedef MultiModel ME;

  typedef System<MG,ME> S;

  typedef ToyFeatures2Multi<ME> F;
  typedef FeaturesInt<F,ME> FI;
  typedef RankAgent<F,ME> RA;

  typedef M2QOptim<S,RA,FI,ME> SPO;

  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen.setType(MLE);
  s.modelEst.setType(MLE);

  int numReps = 1;
  Starts starts(numReps,s.fD.numNodes);
 
  RA ra; // running at the good starting weights

  ra.tp.weights_r.zeros(ra.f.numFeatures);
  ra.tp.weights_r(2) = 1;

  SPO spo;

  R_RA r_ra;

  njm::message("Priority Score: "
  	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts),
  			       ""));

  


//   typedef GravityTimeInfExpCavesModel MG;
  
//   typedef MultiModel ME;

//   typedef System<MG,ME> S;

//   S s;
    
//   int numReps = 96;
//   Starts starts(numReps,s.fD.numNodes);

//   s.reset(starts[0]);
//   s.revert();

//   int i,I=10;
// #pragma omp parallel for			\
//   private(i)					\
//   firstprivate(s)				\
//   shared(I)
//   for(i = 0; i < I; ++i){
//     {
//       System<ME,ME> s2(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.modelEst);
//     }
//   }
    
    
  
  
  
  njm::sett.clean();
  return 0;
}
