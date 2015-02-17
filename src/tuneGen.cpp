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


double TuneGenPA(S & s){
  double trtSize = s.modelGen.tuneTrt(s.fD,s.paramGen);

  s.paramGen_r.trtPre = s.paramGen_r.trtAct = trtSize;
  s.paramEst_r.trtPre = s.paramEst_r.trtAct = trtSize;

  PA pa;
  RP rp;

  return rp.run(s,pa,500,s.fD.finalT);
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  S s;
  s.paramEst_r = s.paramGen_r;
  s.reset();

  njm::message("Tuning Intercept");

  double valNT = TuneGenNT(s);

  njm::message("Tuning Treatment");

  double valPA = TuneGenPA(s);

  njm::message(" intcp: " + njm::toString(s.paramGen_r.intcp,"") +
	       "\n" +
	       "trtPre: " + njm::toString(s.paramGen_r.trtPre,"") +
	       "\n" +
	       "trtAct: " + njm::toString(s.paramGen_r.trtAct,"") +
	       "\n" +
	       " valNT: " + njm::toString(valNT,"") +
	       "\n" +
	       " valPA: " + njm::toString(valPA,""));

  s.paramGen_r.save();
  
  njm::sett.clean();
  
  return 0;
}
