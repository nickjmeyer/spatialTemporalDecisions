#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

//   typedef GravityTimeInfExpCavesModel GM;
//   typedef GravityTimeInfExpCavesParam GP;
//   typedef GM EM;
//   typedef GP EP;

//   typedef System<GM,GP,EM,EP> S;

//   typedef ToyFeatures2<EM,EP> F;
//   typedef RankAgent<F,EM,EP> RA;

//   typedef FeaturesInt<F,EM,EP> FI;
//   typedef M2QOptim<S,RA,FI,EM,EP> OQ;

//   S s;

//   RA ra;
//   OQ oq;

//   s.modelGen.fitType = MLE;
//   s.modelEst.fitType = MLE;

//   s.paramEst_r = s.paramGen_r;
//   s.reset();

//   int i;
//   for(i = 0; i < 1; i++)
//     njm::runif01();
  
//   int t;
//   for(t = 0; t < 10; ++t){
//     if(t >= s.fD.trtStart){
//       ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.paramEst);
//     }
//     s.nextPoint();

//   }

//   std::cout << "value: " << s.value() << std::endl;
  
//   oq.qEval.preCompData(s.sD,s.fD);

//   oq.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst);

//   oq.qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,ra);

//   oq.qEval.buildRD();

// #pragma omp parallel for			\
//   private(i)					\
//   firstprivate(oq)
//   for(i = 0; i < 100; ++i){
//     std::cout << i << std::endl;
//     oq.qEval.solve();
//   }


  std::vector<int> ia,ja;
  std::vector<double> a,b,x;

  njm::fromFile(ia,"/home/nick/research/spatialDecisionMaking/data/toy/grid100/2015-03-25-18-39-25/ia_6.txt");
  njm::fromFile(ja,"/home/nick/research/spatialDecisionMaking/data/toy/grid100/2015-03-25-18-39-25/ja_6.txt");
  njm::fromFile(a,"/home/nick/research/spatialDecisionMaking/data/toy/grid100/2015-03-25-18-39-25/i_6.txt");


  
  njm::sett.clean();
  return 0;
}
