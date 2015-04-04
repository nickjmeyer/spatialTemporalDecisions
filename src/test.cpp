#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel M;
  typedef System<M,M> S;
  // typedef NoTrt<M> N;
  typedef ProximalAgent<M> P;
  // typedef VanillaRunner<S,N> R;

  S s;
  P p;
  // N n;
  // R r;

  // int numReps = 500;
  // Starts starts(numReps,s.fD.numNodes);
  
  // std::cout << r.run(s,n,500,s.fD.finalT,starts) << std::endl;

  Starts starts("startingLocations.txt");

  s.reset(starts[0]);
  s.revert();

  int i;
  for(i = 0; i < 10; ++i){
    if(i >= s.fD.trtStart){
      p.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen);
    }
    s.nextPoint();
  }

  std::cout << "prop inf: " << s.value() << std::endl;

  ToyFeatures2<M> f;
  ToyFeatures2<RangeModel> fr;
  ToyFeatures2<CaveModel> fc;
  ToyFeatures2Multi<MultiModel> fm;

  MultiModel mm;
  RangeModel rm;
  CaveModel cm;

  mm.setType(MLE);
  rm.setType(MLE);
  cm.setType(MLE);
  
  rm.fit(s.sD,s.tD,s.fD,s.dD,0);
  cm.fit(s.sD,s.tD,s.fD,s.dD,0);

  mm.modSel(0);
  mm.getPar()->putPar(s.modelGen.getPar()->getPar());

  mm.modSel(1);
  mm.getPar()->putPar(rm.getPar()->getPar());

  mm.modSel(2);
  mm.getPar()->putPar(cm.getPar()->getPar());
  
  
  f.preCompData(s.sD,s.tD,s.fD,s.dD,s.modelGen);
  fr.preCompData(s.sD,s.tD,s.fD,s.dD,rm);
  fc.preCompData(s.sD,s.tD,s.fD,s.dD,cm);
  fm.preCompData(s.sD,s.tD,s.fD,s.dD,mm);
  
  f.getFeatures(s.sD,s.tD,s.fD,s.dD,s.modelGen);
  fr.getFeatures(s.sD,s.tD,s.fD,s.dD,rm);
  fc.getFeatures(s.sD,s.tD,s.fD,s.dD,cm);
  fm.getFeatures(s.sD,s.tD,s.fD,s.dD,mm);

  std::cout << "===========initial==============" << std::endl;
  std::cout << "InfFeat: " << std::endl;
  std::cout << arma::sum(f.infFeat,0);
  std::cout << arma::sum(fr.infFeat,0);
  std::cout << arma::sum(fc.infFeat,0);
  std::cout << arma::sum(fm.infFeat,0);
  
  std::cout << "NotFeat: " << std::endl;
  std::cout << arma::sum(f.notFeat,0);
  std::cout << arma::sum(fr.notFeat,0);
  std::cout << arma::sum(fc.notFeat,0);
  std::cout << arma::sum(fm.notFeat,0);

  std::cout << std::endl;

  for(i = 0; i < 2; ++i){
    s.tD.p.at(s.sD.notInfec.at(i)) = 1;
    s.tD.a.at(s.sD.infected.at(i)) = 1;
  }

  f.updateFeatures(s.sD,s.tD,s.fD,s.dD,s.modelGen);
  fr.updateFeatures(s.sD,s.tD,s.fD,s.dD,rm);
  fc.updateFeatures(s.sD,s.tD,s.fD,s.dD,cm);
  fm.updateFeatures(s.sD,s.tD,s.fD,s.dD,mm);

  std::cout << "===========updated==============" << std::endl;  
  std::cout << "InfFeat: " << std::endl;
  std::cout << arma::sum(f.infFeat,0);
  std::cout << arma::sum(fr.infFeat,0);
  std::cout << arma::sum(fc.infFeat,0);
  std::cout << arma::sum(fm.infFeat,0);
  
  std::cout << "NotFeat: " << std::endl;
  std::cout << arma::sum(f.notFeat,0);
  std::cout << arma::sum(fr.notFeat,0);
  std::cout << arma::sum(fc.notFeat,0);
  std::cout << arma::sum(fm.notFeat,0);

  std::cout << std::endl;

  for(i = 2; i < 4; ++i){
    s.tD.p.at(s.sD.notInfec.at(i)) = 1;
    s.tD.a.at(s.sD.infected.at(i)) = 1;
  }

  f.updateFeatures(s.sD,s.tD,s.fD,s.dD,s.modelGen);
  fr.updateFeatures(s.sD,s.tD,s.fD,s.dD,rm);
  fc.updateFeatures(s.sD,s.tD,s.fD,s.dD,cm);
  fm.updateFeatures(s.sD,s.tD,s.fD,s.dD,mm);

  std::cout << "===========updated==============" << std::endl;  
  std::cout << "InfFeat: " << std::endl;
  std::cout << arma::sum(f.infFeat,0);
  std::cout << arma::sum(fr.infFeat,0);
  std::cout << arma::sum(fc.infFeat,0);
  std::cout << arma::sum(fm.infFeat,0);
  
  std::cout << "NotFeat: " << std::endl;
  std::cout << arma::sum(f.notFeat,0);
  std::cout << arma::sum(fr.notFeat,0);
  std::cout << arma::sum(fc.notFeat,0);
  std::cout << arma::sum(fm.notFeat,0);

  std::cout << std::endl;

  for(i = 4; i < 6; ++i){
    s.tD.p.at(s.sD.notInfec.at(i)) = 1;
    s.tD.a.at(s.sD.infected.at(i)) = 1;
  }

  f.updateFeatures(s.sD,s.tD,s.fD,s.dD,s.modelGen);
  fr.updateFeatures(s.sD,s.tD,s.fD,s.dD,rm);
  fc.updateFeatures(s.sD,s.tD,s.fD,s.dD,cm);
  fm.updateFeatures(s.sD,s.tD,s.fD,s.dD,mm);

  std::cout << "===========updated==============" << std::endl;  
  std::cout << "InfFeat: " << std::endl;
  std::cout << arma::sum(f.infFeat,0);
  std::cout << arma::sum(fr.infFeat,0);
  std::cout << arma::sum(fc.infFeat,0);
  std::cout << arma::sum(fm.infFeat,0);
  
  std::cout << "NotFeat: " << std::endl;
  std::cout << arma::sum(f.notFeat,0);
  std::cout << arma::sum(fr.notFeat,0);
  std::cout << arma::sum(fc.notFeat,0);
  std::cout << arma::sum(fm.notFeat,0);
  
  
  njm::sett.clean();
  return 0;
}
