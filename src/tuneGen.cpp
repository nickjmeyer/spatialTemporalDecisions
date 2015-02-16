#include "tuneGen.hpp"

void TuneGenNT(){
  S s;
  NT nt;
  RN rn;

  double goal = 0.7;
  int numReps = 500;
  int numYears = s.fD.finalT;
  double tol = 0.001;

  double par = s.paramGen_r.intcp;
  double val = rn.run(s,nt,numReps,numYears);
  double add = 1.0, scale = .975;
  int above = int(val > goal);
  int iter = 0;

  printf("Iter: %05d  >>>  Current value: %012.6f  @  %08.4f\r",
	 ++iter, val, par);

  while(std::abs(val - goal) > tol){
    if(val > goal){
      if(!above)
	add*=scale;
      
      par -= add;
      s.paramGen_r.intcp = par;
      s.paramEst_r.intcp = par;

      above = 1;
    }
    else{
      if(above)
	add*=scale;

      par += add;
      s.paramGen_r.intcp = par;
      s.paramEst_r.intcp = par;
      
      above = 0;
    }

    val = rn.run(s,nt,numReps,numYears);
    printf("Iter: %05d  >>>  Current value: %012.6f  @  %08.4f\r",
	   ++iter, val, par);
    fflush(stdout);
  }
}


void TuneGenMA(){
  S s;
  s.paramEst_r = s.paramGen_r;
  MA ma;
  RM rm;

  double goal = 0.5;
  int numReps = 500;
  int numYears = s.fD.finalT;
  double tol = 0.001;

  double par = s.paramGen_r.trtPre;
  double val = rm.run(s,ma,numReps,numYears);
  double add = 1.0, scale = .975;
  int above = int(val > goal);
  int iter = 0;

  printf("Iter: %05d  >>>  Current value: %012.6f  @  %08.4f\r",
	 ++iter, val, par);

  while(std::abs(val - goal) > tol){
    if(val > goal){
      if(!above)
	add*=scale;
      
      par += add;
      s.paramGen_r.trtPre = par;
      s.paramEst_r.trtPre = par;
      s.paramGen_r.trtAct = par;
      s.paramEst_r.trtAct = par;

      above = 1;
    }
    else{
      if(above)
	add*=scale;

      par -= add;
      s.paramGen_r.trtPre = par;
      s.paramEst_r.trtPre = par;
      s.paramGen_r.trtAct = par;
      s.paramEst_r.trtAct = par;
      
      above = 0;
    }

    val = rm.run(s,ma,numReps,numYears);
    printf("Iter: %05d  >>>  Current value: %012.6f  @  %08.4f\r",
	   ++iter, val, par);
    fflush(stdout);
  }
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  njm::message("Tuning Intercept");

  TuneGenNT();

  njm::message("Tuning Treatment");

  TuneGenMA();

  
  njm::sett.clean();
  
  return 0;
}
