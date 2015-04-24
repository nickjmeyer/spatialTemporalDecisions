#include "modelUnits.hpp"


void test(const std::string & name, const int cond){
  if(cond)
    printf("%32s: passed\n",name.c_str());
  else
    printf("%32s: failed ***\n",name.c_str());
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  // typedef ModelTimeExpCaves MG;
  typedef ModelTime MG;
  // typedef ModelGravity MG;
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
  std::vector<double> par = s.modelGen.getPar();
  std::cout << "par: " << njm::toString(par," ","\n");
  s.modelEst_r.putPar(par.begin());
  s.modelEst_r.setFill(s.sD,s.tD,s.fD,s.dD);
  s.revert();

  std::vector<double> infProbs;
  // s.modelGen.setFill(s.sD,s.tD,s.fD,s.dD);
  // s.modelGen.setQuick(s.sD,s.tD,s.fD,s.dD);
  s.modelGen.infProbs(s.sD,s.tD,s.fD,s.dD);
  infProbs = s.modelGen.infProbs();
  std::cout << njm::toString(s.sD.infected," ","\n");
  std::cout << njm::toString(s.modelGen.getPar()," ","\n");
  std::cout << std::accumulate(infProbs.begin(),infProbs.end(),0.0)
	    << std::endl
	    << s.modelGen.oneOnOne(0,89,s.fD.numNodes)
	    << std::endl
	    << s.modelGen.oneOnOne(1,89,s.fD.numNodes)
	    << std::endl
	    << s.modelGen.oneOnOne(2,89,s.fD.numNodes)
	    << std::endl;

  s.nextPoint();

  s.modelGen.infProbs(s.sD,s.tD,s.fD,s.dD);
  infProbs = s.modelGen.infProbs();
  std::cout << njm::toString(s.sD.infected," ","\n");
  std::cout << njm::toString(s.modelGen.getPar()," ","\n");
  std::cout << std::accumulate(infProbs.begin(),infProbs.end(),0.0)
	    << std::endl
	    << s.modelGen.oneOnOne(0,89,s.fD.numNodes)
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
