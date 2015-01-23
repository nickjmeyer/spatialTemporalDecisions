#include "trainLambda.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam_r = s.genParam_r;
  s.reset();

  OptimRunner<System,RankToyAgent,GravityModel,GravityParam,M2NmOptim> pR;
  RankToyAgent<GravityModel,GravityParam> rA;
  M2NmOptim<System,RankToyAgent,GravityModel,GravityParam> qO;

  int i,minLambda=0;
  double val,minVal=1.0;
  for(i=5000; i<15001; i+=1000){
    qO.qEval.tp.lambda=i;
    val = pR.run(s,rA,qO,150,s.fD.finalT);
    njm::message("lambda: " + njm::toString(i," ",6,0)
		 + "  -->  " + njm::toString(val,"\n",6,4));
    if(val < minVal){
      minLambda = i;
      minVal = val;
    }
  }

  for(i=minLambda-2000; i<minLambda+2001; i+=100){
    qO.qEval.tp.lambda=i;
    val = pR.run(s,rA,qO,150,s.fD.finalT);
    njm::message("lambda: " + njm::toString(i," ",6,0)
		 + "  -->  " + njm::toString(val,"\n",6,4));
    if(val < minVal){
      minLambda = i;
      minVal = val;
    }
  }
						  
  njm::message("min :: lambda: " + njm::toString(minLambda," ",6,0)
	       + "  -->  " + njm::toString(minVal,"\n",6,4));
  

  
  
  njm::sett.clean();
  return 0;
}
