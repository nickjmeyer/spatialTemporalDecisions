#include "tuneGen.hpp"

double TuneGenNT(S & s){
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

  printf("Iter: %05d  >>>  Current value: %08.6f  @  %08.4f\r",
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
    printf("Iter: %05d  >>>  Current value: %08.6f  @  %08.4f\r",
	   ++iter, val, par);
    fflush(stdout);
  }

  return(val);
}


double TuneGenMA(S & s){
  // MA ma;
  // RM rm;
  RankToyAgent<ToyFeatures2<EM,EP>,EM,EP> ma;
  VanillaRunnerNS<S,RankToyAgent<ToyFeatures2<EM,EP>,EM,EP> >rm;

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
      s.reset();

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
      s.reset();
      
      above = 0;
    }

    val = rm.run(s,ma,numReps,numYears);
    printf("Iter: %05d  >>>  Current value: %012.6f  @  %08.4f\r",
	   ++iter, val, par);
    fflush(stdout);
  }
  return(val);
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  S s;
  s.paramEst_r = s.paramGen_r;
  s.reset();

  njm::message("Tuning Intercept");

  double valNT = TuneGenNT(s);

  njm::message("Tuning Treatment");

  double valMA = TuneGenMA(s);

  njm::message(" intcp: " + njm::toString(s.paramGen_r.intcp,"") +
	       "\n" +
	       "trtPre: " + njm::toString(s.paramGen_r.trtPre,"") +
	       "\n" +
	       "trtAct: " + njm::toString(s.paramGen_r.trtAct,"") +
	       "\n" +
	       " valNT: " + njm::toString(valNT,"") +
	       "\n" +
	       " valMA: " + njm::toString(valMA,""));

  // s.paramGen_r.save();
  
  njm::sett.clean();
  
  return 0;
}
