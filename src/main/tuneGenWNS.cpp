#include <gflags/gflags.h>
#include "tuneGenWNS.hpp"

DEFINE_string(srcDir,"","Path to source directory");

template <class S, class MA, class RM, class NT, class RN>
double TuneGenMA(S & s, const int numReps, const Starts & starts){
  NT nt;
  RN rn;

  MA ma;
  RM rm;

  int numYears = s.fD.finalT;

  double atTrtStart = rn.run(s,nt,numReps,s.fD.trtStart,starts).smean();
  double atFinalT = rn.run(s,nt,numReps,numYears,starts).smean();

  double goal = atTrtStart + 0.05*(atFinalT - atTrtStart);
  njm::message("Goal: " + njm::toString(goal,""));
  double tol = 1e-3;

  std::vector<double> par;
  double trt = 1.0;

  s.modelGen_r.setPar(std::vector<std::string>({"trtAct","trtPre"}),trt);
  s.modelGen_r.save();
  s = S("obsData.txt");

  double val = rm.run(s,ma,numReps,numYears,starts).smean();
  double scale = 1.1, shrink = .9;
  int above = int(val > goal);
  int iter = 0;


  printf("Iter: %05d  >>>  Current value: %08.6f  ===  Current Trt: %08.6f\r",
	 ++iter, val, trt);

  while(std::abs(val - goal) > tol){
    if(val > goal){
      if(!above)
	scale*=shrink;

      trt *= 1.0 + scale;

      above = 1;
    }
    else{
      if(above)
	scale*=shrink;

      trt *= 1.0/(1.0 + scale);

      above = 0;
    }


    s.modelGen_r.setPar(std::vector<std::string>({"trtAct","trtPre"}),trt);
    s.modelGen_r.save();
    s = S("obsData.txt");

    // std::cout << "par: " << njm::toString(par," ","\n");
    // par = s.modelGen.getPar({"trtAct","trtPre"});
    // std::cout << "par: " << njm::toString(par," ","\n");


    val = rm.run(s,ma,numReps,numYears,starts).smean();
    printf("Iter: %05d  >>>  Current value: %08.6f  ===  Current Trt: %08.6f\r",
	   ++iter, val, trt);
    fflush(stdout);
  }

  s.modelGen_r.save();

  printf("\n");

  njm::message("Est. goal: " + njm::toString(val,""));

  njm::message("par: " + njm::toString(s.modelGen_r.getPar()," ",""));

  return(val);
}


// template <class S, class PA, class RP>
// double TuneGenPA(S & s,const int numReps, const Starts & starts){
//   double trtSize = s.modelGen.tuneTrt(s.fD);

//   putActTrt(trtSize,s.modelGen_r,s.fD);
//   putPreTrt(trtSize,s.modelGen_r,s.fD);
//   putActTrt(trtSize,s.modelEst_r,s.fD);
//   putPreTrt(trtSize,s.modelEst_r,s.fD);

//   PA pa;
//   RP rp;

//   return rp.run(s,pa,numReps,s.fD.finalT,starts).smean();
// }


int main(int argc, char ** argv){
  gflags::ParseCommandLineFlags(&argc,&argv,true);
  njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

  {
    // typedef ModelTimeExpCavesGDistTrendPowCon MG;
    typedef Model2GravityGDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;
    typedef NoTrt<ME> NT;
    typedef ProximalGDistAgent<ME> PA;
    typedef MyopicAgent<ME> MA;

    typedef AllAgent<ME> AA;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef VanillaRunnerNS<S,NT> RN;
    typedef VanillaRunnerNS<S,PA> RP;
    typedef VanillaRunnerNS<S,MA> RM;
    typedef VanillaRunnerNS<S,RA> RR;

    typedef VanillaRunnerNS<S,AA> R_AA;

    S s("obsData.txt");
    s.modelEst_r = s.modelGen_r;
    s.revert();

    int numReps = 500;
    Starts starts("startingLocations.txt");

    NT nt;
    MA ma;
    PA pa;
    RP rp;

    RN rn;
    RA ra;
    RM rm;
    RR rr;
    // ra.reset();

    double valNT = rn.run(s,nt,numReps,s.fD.finalT,starts).smean();

    njm::message("Tuning Treatment");

    double valAA = TuneGenMA<S,AA,R_AA,NT,RN>(s,numReps,starts);

    double valMA = rm.run(s,ma,numReps,s.fD.finalT,starts).smean();

    double valPA = rp.run(s,pa,numReps,s.fD.finalT,starts).smean();

    double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).smean();

    njm::message(" valNT: " + njm::toString(valNT,"") +
		 "\n" +
		 " valPA: " + njm::toString(valPA,"") +
		 "\n" +
		 " valMA: " + njm::toString(valMA,"") +
		 "\n" +
		 " valRA: " + njm::toString(valRA,"") +
		 "\n" +
		 " valAA: " + njm::toString(valAA,""));


    std::vector<double> par = s.modelGen_r.getPar();

    double priorMeanTrt = (s.modelGen_r.getPar({"trtAct"})[0]
			   + s.modelGen_r.getPar({"trtPre"})[0])/2.0;
    priorMeanTrt *= 4.0;

    // write prior mean of treatment effect
    njm::toFile(priorMeanTrt,njm::sett.srcExt("priorTrtMean.txt"),
		std::ios_base::out);
  }

  njm::sett.clean();

  return 0;
}
