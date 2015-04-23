#include "test.hpp"
#include <omp.h>

void fn0(const SimData & sD){
  std::cout << njm::toString(sD.notInfec," ","\n");
}

void fn1(const TrtData & tD){
  std::cout << njm::toString(tD.a," ","\n");
}

void fn2(const FixedData & fD){
  std::cout << njm::toString(fD.numCovar,"\n");
}

void fn3(const DynamicData & dD){
  std::cout << njm::toString("blah","\n");
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelGravity MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  typedef ToyFeatures2<ME> F;
  typedef RankAgent<F,ME> RA;
  typedef PlainRunner<S,RA> PR;
  typedef NoTrt<ME> NT;

  S s;
  RA ra;
  PR pr;
  NT nt;

  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  std::vector<double> par = s.modelGen_r.getPar();
  s.modelEst_r.putPar(par.begin());

  int numReps = 100;
  Starts starts("startingLocations.txt");

  s.reset(starts[0]);
  s.revert();

  njm::timer.start("everything");
  njm::message(pr.run(s,ra,numReps,s.fD.finalT));
  njm::timer.stop("everything");
  // // std::cout << "\n\nFuckers\n\n";
  // // njm::message(pr.run(s,ra,1,s.fD.finalT));

  

  njm::sett.clean();
  return 0;
}
