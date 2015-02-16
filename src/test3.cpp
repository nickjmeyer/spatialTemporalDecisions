#include "test3.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef GravityModel EM;
  typedef GravityParam EP;

  typedef System<GM,GP,EM,EP> S;

  typedef ToyFeatures2<EM,EP> F;

  typedef RankToyAgent<F,EM,EP> RA;

  typedef PlainRunner<S,RA> PR;

  typedef M1SgdOptim<S,RA,EM,EP> M1;

  typedef OptimRunner<S,RA,M1> OR;
  
  S s;
  RA ra;
  PR pr;
  M1 m1;
  OR o1;
  
  s.paramEst_r = s.paramGen_r;
  s.reset();

  njm::message(pr.run(s,ra,600,15));

  std::vector<double> aVals;
  aVals.push_back(15);
  aVals.push_back(30);
  aVals.push_back(50);
  aVals.push_back(100);
  aVals.push_back(150);
  aVals.push_back(200);
  std::vector<double>::const_iterator aIt;

  std::vector<int> mcRepsVals;
  mcRepsVals.push_back(25);
  mcRepsVals.push_back(50);
  mcRepsVals.push_back(100);
  std::vector<int>::const_iterator mcRepsIt;

  double value;
  for(mcRepsIt = mcRepsVals.begin(); mcRepsIt != mcRepsVals.end(); mcRepsIt++){
    for(aIt = aVals.begin(); aIt != aVals.end(); aIt++){
      m1.tp.mcReps = (*mcRepsIt);
      m1.tp.a = (*aIt);
      
      value = o1.run(s,ra,m1,600,15);
    
      printf("reps: %3d, a: %3.0f  >>>  %6.4f\n",*mcRepsIt,*aIt,value);
    }
  }

  njm::sett.clean();
  return 0;
}
