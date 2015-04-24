#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCaves MG;
  typedef ModelTime ME;
  typedef System<MG,ME> S;
  typedef ToyFeatures2<ME> F;
  typedef RankAgent<F,ME> RA;
  // typedef PlainRunner<S,RA> PR;
  typedef FitOnlyRunner<S,RA> FR;

  S s;
  RA ra;
  // PR pr;
  FR fr;

  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  // std::vector<double> par = s.modelGen_r.getPar();
  // s.modelEst_r.putPar(par.begin());

  int numReps = 20;
  Starts starts("startingLocations.txt");

  omp_set_num_threads(1);

  s.reset(starts[0]);
  s.revert();

  njm::timer.start("everything");
  // njm::message(pr.run(s,ra,numReps,s.fD.finalT));
  njm::message(fr.run(s,ra,numReps,s.fD.finalT,starts));
  njm::timer.stop("everything");

  njm::sett.clean();
  return 0;
}
