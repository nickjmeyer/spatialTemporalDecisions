#include "test2.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  typedef GravityModel MG;
  typedef GravityParam PG;
  
  typedef GravityModel ME;
  typedef GravityParam PE;

  typedef System<MG,PG,ME,PE> S;
  
  typedef ToyFeatures2<ME,PE> F;
  
  typedef RankToyAgent<F,ME,PE> AR;

  typedef VanillaRunner<S,AR> R_AR;

  // system
  S s;
  s.modelEst = s.modelGen;
  s.paramEst_r = s.paramGen_r;
  s.paramEst = s.paramGen;

  // agents
  AR ar;

  // runners
  R_AR r_ar;


  std::vector<double> jitterVals;
  jitterVals.push_back(0.00);
  jitterVals.push_back(0.05);
  jitterVals.push_back(0.10);
  jitterVals.push_back(0.25);
  jitterVals.push_back(0.50);
  jitterVals.push_back(0.75);
  jitterVals.push_back(1.00);
  jitterVals.push_back(1.25);
  jitterVals.push_back(1.50);
  jitterVals.push_back(1.75);
  jitterVals.push_back(2.00);
  jitterVals.push_back(2.25);
  jitterVals.push_back(2.50);
  jitterVals.push_back(2.75);
  jitterVals.push_back(3.00);
  jitterVals.push_back(4.00);
  jitterVals.push_back(4.25);
  jitterVals.push_back(4.50);
  jitterVals.push_back(4.75);
  jitterVals.push_back(5.00);


  std::vector<double> values;

  int mcReps=2000,numPoints = s.fD.finalT;

  for(int i=0; i<(int)jitterVals.size(); i++){
    printf("done %03d / %03d\r",i,(int)jitterVals.size());
    fflush(stdout);
    
    ar.tp.jitter = jitterVals.at(i);
    
    values.push_back(r_ar.run(s,ar,mcReps,numPoints));
  }

  for(int i=0; i<(int)jitterVals.size(); i++){
    printf("jitter: %05.3f  -->  value: %06.4f\n",
	   jitterVals.at(i),values.at(i));
  }

  
  njm::sett.clean();
  return 0;
}
