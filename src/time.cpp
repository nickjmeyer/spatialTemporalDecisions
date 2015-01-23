#include "time.hpp"
#include "omp.h"

int main(){

  System<GravityModel,GravityParam> s;
  
  RankAgent<GravityModel,GravityParam> rA;
  PlainRunner<System,RankAgent,GravityModel,GravityParam> rR;
  
  // RankAgent<GravityModel,GravityParam> rA;
  // PlainRunner<System<GravityModel,GravityParam>,
  // 	      RankAgent<GravityModel,GravityParam> > rR;

  // NoTrt<GravityModel,GravityParam> rA;
  // PlainRunner<System<GravityModel,GravityParam>,
  // 	      NoTrt<GravityModel,GravityParam> > rR;

  
  s.estParam_r = s.genParam_r;
  s.reset();

  
  
  int threads=omp_get_max_threads();
  int i,N=(threads < 32 ? 1 : threads),reps=20;
  njm::message("threads: " + njm::toString(threads));
  njm::message("N: " + njm::toString(N));

  int tick = std::time(NULL);
  if(N>1){
#pragma omp parallel for num_threads(threads)	\
  private(i)					\
  firstprivate(rR,s,rA,reps)			\
  shared(N)
    for(i=0; i<N; i++)
    njm::message(rR.run(s,rA,reps,s.fD.finalT));
  }
  else
    njm::message(rR.run(s,rA,reps,s.fD.finalT));    
  int tock = std::time(NULL);
    
  double elapsed = tock - tick;
  
  njm::message("Time elapsed: " + njm::toString(elapsed,"") + " seconds\n");
  njm::message("Time elapsed: "
  	       + njm::toString(elapsed/((double) N),"")
  	       + " seconds/thread\n");
  njm::message("Time elapsed: "
  	       + njm::toString(elapsed/((double) reps),"")
  	       + " seconds/run\n");
  njm::message("Time elapsed: "
  	       + njm::toString(elapsed/((double) reps*N),"")
  	       + " seconds/(thread*run)\n");
  
  elapsed/=3600.0;
  
  njm::message("Time elapsed: " + njm::toString(elapsed,"") + " hours\n");
  njm::message("Time elapsed: "
  	       + njm::toString(elapsed/((double) N),"")
  	       + " hours/thread\n");
  njm::message("Time elapsed: "
  	       + njm::toString(elapsed/((double) reps),"")
  	       + " hours/run\n");
  njm::message("Time elapsed: "
  	       + njm::toString(elapsed/((double) N*reps),"")
  	       + " hours/(thread*run)\n");
  
  return 0;
}
