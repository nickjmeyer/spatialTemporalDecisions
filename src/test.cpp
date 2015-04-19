#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);



  // typedef GravityModel MG;
  // typedef MG ME;
  // typedef System<MG,ME> S;
  // typedef ToyFeatures2<MG> F;
  // typedef RankAgent<F,ME> RA;
  // typedef PlainRunner<S,RA> PR;

  // S s;
  // RA ra;
  // PR pr;

  // s.modelGen_r.setType(MLE);
  // s.modelEst_r.setType(MLE);

  // s.modelEst_r = s.modelGen_r;
  

  // int numReps = 100;
  // Starts starts("startingLocations.txt");

  // s.reset(starts[0]);

  // njm::timer.start("everything");
  // njm::message(pr.run(s,ra,numReps,s.fD.finalT));
  // njm::timer.stop("everything");
  // // std::cout << "\n\nFUCKERS\n\n";
  // // njm::message(pr.run(s,ra,1,s.fD.finalT));



  ////////////////////////////////////////

  typedef GravityTimeInfExpCavesModel MG;
  // typedef GravityModel MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  // typedef NoTrt<MG> NT;
  // typedef PlainRunner<S,NT> PR;
  typedef ToyFeatures2<MG> F;
  typedef RankAgent<F,ME> RA;
  // typedef PlainRunner<S,RA> PR;
  typedef FitOnlyRunner<S,RA> FR;

  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  // NT nt;
  RA ra;
  // PR pr;
  FR pr;
  
  int numReps = 1;
  // omp_set_num_threads(2);
  Starts starts("startingLocations.txt");

  s.reset(starts[0]);
  std::vector<double> par = s.modelGen.getPar()->getPar();
  std::cout << "par: " << njm::toString(par," ","\n");
  s.modelEst_r.getPar()->putPar(par);
  s.revert();

  std::vector<double> infProbs;
  s.modelGen.infProbs(s.sD,s.tD,s.fD,s.dD);
  infProbs = s.modelGen.getPar()->getInfProbs();
  std::cout << njm::toString(s.sD.infected," ","\n");
  std::cout << njm::toString(s.modelGen.getPar()->getPar()," ","\n");
  std::cout << std::accumulate(infProbs.begin(),infProbs.end(),0.0)
	    << std::endl
	    << s.modelGen.oneOnOne(0,89,s.sD,s.tD,s.fD,s.dD)
	    << std::endl
	    << s.modelGen.oneOnOne(1,89,s.sD,s.tD,s.fD,s.dD)
	    << std::endl
	    << s.modelGen.oneOnOne(2,89,s.sD,s.tD,s.fD,s.dD)
	    << std::endl;

  s.nextPoint();

  s.modelGen.infProbs(s.sD,s.tD,s.fD,s.dD);
  infProbs = s.modelGen.getPar()->getInfProbs();
  std::cout << njm::toString(s.sD.infected," ","\n");
  std::cout << njm::toString(s.modelGen.getPar()->getPar()," ","\n");
  std::cout << std::accumulate(infProbs.begin(),infProbs.end(),0.0)
	    << std::endl
	    << s.modelGen.oneOnOne(0,89,s.sD,s.tD,s.fD,s.dD)
	    << std::endl;

  F f;
  f.preCompData(s.sD,s.tD,s.fD,s.dD,s.modelEst);
  f.getFeatures(s.sD,s.tD,s.fD,s.dD,s.modelEst);

  std::cout << arma::sum(f.infFeat,0)
	    << arma::sum(f.notFeat,0);

  s.tD.a.at(s.sD.infected.at(0)) = 1;
  s.tD.a.at(s.sD.infected.at(1)) = 1;

  s.tD.p.at(s.sD.notInfec.at(0)) = 1;
  s.tD.p.at(s.sD.notInfec.at(1)) = 1;

  f.updateFeatures(s.sD,s.tD,s.fD,s.dD,s.modelEst);

  std::cout << arma::sum(f.infFeat,0)
	    << arma::sum(f.notFeat,0);

  // std::cout << "run: " << pr.run(s,nt,numReps,s.fD.finalT)
  // 	    << std::endl;
  std::cout << "run: " << pr.run(s,ra,numReps,s.fD.finalT,starts)
  	    << std::endl;

  

  njm::sett.clean();
  return 0;
}
