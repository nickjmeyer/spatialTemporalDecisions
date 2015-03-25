#include "test.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel GM;
  typedef GravityTimeInfExpCavesParam GP;
  typedef GM EM;
  typedef GP EP;

  typedef System<GM,GP,EM,EP> S;

  typedef ToyFeatures2<EM,EP> F;
  typedef RankAgent<F,EM,EP> RA;

  typedef FeaturesInt<F,EM,EP> FI;
  typedef M2QOptim<S,RA,FI,EM,EP> OQ;

  S s;

  RA ra;
  OQ oq;

  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  s.paramEst_r = s.paramGen_r;
  s.reset();

  int i;
  for(i = 0; i < 0; i++)
    njm::runif01();
  
  int t;
  for(t = 0; t < s.fD.finalT; ++t){
    if(t >= s.fD.trtStart){
      ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.paramEst);
    }
    s.nextPoint();

  }

  std::cout << "value: " << s.value() << std::endl;
  
  oq.qEval.preCompData(s.sD,s.fD);

  oq.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst);

  oq.qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,ra);

  oq.qEval.buildRD();

  njm::timer.start("SparseLU");
  for(i = 0; i < 10; ++i){
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.compute(oq.qEval.DtD + oq.qEval.tp.lambda*oq.qEval.P);
    oq.qEval.beta = solver.solve(oq.qEval.mDtR);
  }
  njm::timer.stop("SparseLU");

  njm::timer.start("SuperLU");
  for(i = 0; i < 10; ++i){
    Eigen::SuperLU<Eigen::SparseMatrix<double> > solver;
    solver.compute(oq.qEval.DtD + oq.qEval.tp.lambda*oq.qEval.P);
    oq.qEval.beta = solver.solve(oq.qEval.mDtR);
  }
  njm::timer.stop("SuperLU");

  // njm::timer.start("PardisoLU");
  // for(i = 0; i < 10; ++i){
  //   Eigen::PardisoLU<Eigen::SparseMatrix<double> > solver;
  //   // solver.compute(oq.qEval.DtD + oq.qEval.tp.lambda*oq.qEval.P);
  //   // oq.qEval.beta = solver.solve(oq.qEval.mDtR);
  // }
  // njm::timer.stop("PardisoLU");
  
  njm::sett.clean();
  return 0;
}
