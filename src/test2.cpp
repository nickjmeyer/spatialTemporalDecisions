#include "test2.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  typedef GravityModel MG;
  typedef GravityParam PG;
  
  typedef MG ME;
  typedef PG PE;

  typedef System<MG,PG,ME,PE> S;
  
  typedef ToyFeatures2<ME,PE> F;
  
  typedef RankToyAgent<F,ME,PE> AR;

  typedef M1SpOptim<S,AR,ME,PE> SPO;
  typedef M1SgdOptim<S,AR,ME,PE> SGDO;

  typedef OptimRunner<S,AR,SPO> SPR;
  typedef OptimRunner<S,AR,SGDO> SGDR;

  // system
  S s;
  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;
  s.reset();

  AR ar;

  SPO spo;
  spo.tp.tune = 0;
  SGDO sgdo;

  SPR spr;
  SGDR sgdr;

  int numYears = 15;
  int numReps = 300;

  // omp_set_num_threads(1);
  
  njm::message(spr.run(s,ar,spo,numReps,numYears));
  njm::message(sgdr.run(s,ar,sgdo,numReps,numYears));


  // s.paramEst_r = s.paramGen_r;

  // PlainRunner<S,AR> pr;

  // pr.run(s,ar,66,15);
  

  
  njm::sett.clean();
  return 0;
}
