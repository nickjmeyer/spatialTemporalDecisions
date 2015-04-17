#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);



  typedef GravityTimeInfExpCavesModel MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  typedef ToyFeatures2<MG> F;
  typedef RankAgent<F,ME> RA;
  typedef PlainRunner<S,RA> PR;

  S s;
  RA ra;
  PR pr;


  // int numReps = 1000;
  // Starts starts("startingLocations.txt");

  // s.reset(starts[0]);

  // pr.run(s,ra,numReps,s.fD.finalT);

  ParamBase * pB = new ParamIntercept(s.fD);

  std::cout << njm::toString(pB->getPar()," ","\n");

  std::vector<double> pars = {1};
  std::vector<double>::iterator beg = pars.begin();

  pB->putPar(beg);

  std::cout << njm::toString(pB->getPar()," ","\n");

  std::cout << (beg == pars.end()) << std::endl;


  std::vector<double> probs = {0,1,2};

  std::cout << njm::toString(probs," ","\n");

  pB->setFill(probs);
  
  std::cout << njm::toString(probs," ","\n");

  pars = {2};
  pB->putPar(pars.begin());
  pB->updateFill(probs);
  
  std::cout << njm::toString(probs," ","\n");

  probs = {0,1,2};

  pars = {1};
  pB->putPar(pars.begin());
  pB->setFill(probs);
  std::cout << njm::toString(probs," ","\n");

  pars = {2};
  pB->putPar(pars.begin());
  pB->setFill(probs);
  std::cout << njm::toString(probs," ","\n");

  pB->updateFill(probs);
  std::cout << njm::toString(probs," ","\n");

  pars = {-1};
  pB->putPar(pars.begin());
  pB->updateFill(probs);
  std::cout << njm::toString(probs," ","\n");
  
  

  
  
  
  delete pB;


  njm::sett.clean();
  return 0;
}
