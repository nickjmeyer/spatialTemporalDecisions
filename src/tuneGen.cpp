#include "tuneGen.hpp"

template <class S, class NT,class RN>
double TuneGenNT(S & s){
  NT nt;
  RN rn;

  double goal = 0.7;
  int numReps = 500;
  int numYears = s.fD.finalT;
  double tol = 0.001;

  std::vector<double> par = s.paramGen_r.getPar();
  double power = s.paramGen_r.power;
  double val = rn.run(s,nt,numReps,numYears);
  double scale = 3.0, shrink = .975;
  int above = int(val > goal);
  int iter = 0;

  printf("Iter: %05d  >>>  Current value: %08.6f\r",
	 ++iter, val);

  while(std::abs(val - goal) > tol){
    if(val > goal){
      if(!above)
	scale*=shrink;
      
      std::for_each(par.begin(),par.end(),
		    [&scale](double & x){x*= 1.0 + scale;});
      
      s.paramGen_r.putPar(par);
      s.paramGen_r.power = power;
      s.paramEst_r.putPar(par);
      s.paramEst_r.power = power;
      
      s.reset();

      above = 1;
    }
    else{
      if(above)
	scale*=shrink;

      std::for_each(par.begin(),par.end(),
		    [&scale](double & x){x*= 1.0/(1.0 + scale);});
      
      s.paramGen_r.putPar(par);
      s.paramGen_r.power = power;
      s.paramEst_r.putPar(par);
      s.paramEst_r.power = power;
      
      s.reset();
      
      above = 0;
    }

    val = rn.run(s,nt,numReps,numYears);
    printf("Iter: %05d  >>>  Current value: %08.6f\r", ++iter, val);
    fflush(stdout);
  }

  return(val);
}


template <class S, class PA, class RP>
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

  {
    typedef GravityModel GM;
    typedef GravityParam GP;
    typedef GM EM;
    typedef GP EP;

    typedef System<GM,GP,EM,EP> S;
    typedef NoTrt<EM,EP> NT;
    typedef ProximalAgent<EM,EP> PA;

    typedef VanillaRunnerNS<S,NT> RN;
    typedef VanillaRunnerNS<S,PA> RP;

    S s;
    s.paramEst_r = s.paramGen_r;
    s.reset();

    njm::message("Tuning Intercept");

    double valNT = TuneGenNT<S,NT,RN>(s);

    njm::message("Tuning Treatment");

    double valPA = TuneGenPA<S,PA,RP>(s);

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
  }

  
  {
    typedef GravityTimeInfModel GM;
    typedef GravityTimeInfParam GP;
    typedef GM EM;
    typedef GP EP;

    typedef System<GM,GP,EM,EP> S;
    typedef NoTrt<EM,EP> NT;
    typedef ProximalAgent<EM,EP> PA;

    typedef VanillaRunnerNS<S,NT> RN;
    typedef VanillaRunnerNS<S,PA> RP;

    S s;
    s.paramEst_r = s.paramGen_r;
    s.reset();

    njm::message("Tuning Intercept");

    double valNT = TuneGenNT<S,NT,RN>(s);

    njm::message("Tuning Treatment");

    double valPA = TuneGenPA<S,PA,RP>(s);

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


    double priorMeanTrt = (s.paramGen_r.trtPre + s.paramGen_r.trtAct)/2.0;
    priorMeanTrt *= 4.0;

    // write prior mean of treatment effect
    njm::toFile(priorMeanTrt,njm::sett.srcExt("priorTrtMean.txt"),
		std::ios_base::out);
  }


  
  njm::sett.clean();
  
  return 0;
}
