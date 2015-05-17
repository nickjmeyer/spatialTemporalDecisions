#include "tuneGen.hpp"

// template <class M>
// double getAlpha(const M & m,
// 		const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   return par.at(1 + fD.numCovar);
// }

// template <class M>
// void putAlpha(const double & alpha,
// 	      M & m,
// 	      const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   par.at(1 + fD.numCovar) = alpha;
//   m.putPar(par.begin());
// }

// template <class M>
// double getPower(const M & m,
// 		const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   return par.at(1 + fD.numCovar + 1);
// }

// template <class M>
// void putPower(const double & power,
// 	      M & m,
// 	      const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   par.at(1 + fD.numCovar + 1) = power;
//   m.putPar(par.begin());
// }

// template <class M>
// double getActTrt(const M & m,
// 		 const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   return par.at(par.size() - 2);
// }

// template <class M>
// void putActTrt(const double & trt,
// 	       M & m,
// 	       const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   par.at(par.size() - 2) = trt;
//   m.putPar(par.begin());
// }


// template <class M>
// double getPreTrt(const M & m,
// 		 const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   return par.at(par.size() - 1);
// }

// template <class M>
// void putPreTrt(const double & trt,
// 	       M & m,
// 	       const FixedData & fD){
//   std::vector<double> par = m.getPar();
//   par.at(par.size() - 1) = trt;
//   m.putPar(par.begin());
// }



double getDPow(const double & power, const double & alpha,
	       const std::vector<double> & caves){
  double meanCaves = std::accumulate(caves.begin(),caves.end(),0);
  meanCaves /= double(caves.size());

  double dPow = std::log(2.0)*std::pow(meanCaves,2.0*power)/alpha + 1.0;
  dPow = std::log(dPow);
  dPow /= std::log(2.0);

  return(dPow);
  // return(1.0);
}


void rescaleD(const double & pastScale, const double & currScale,
	      std::vector<double> & d){
  double scale = currScale/pastScale;
  std::for_each(d.begin(),d.end(),[&scale](double & x){x=std::pow(x,scale);});
}



template <class S, class NT,class RN, class MG>
double TuneGenNT(S & s, const int numReps, const Starts & starts){
  NT nt;
  RN rn;

  double goal = 0.7;
  njm::message("Goal: " + njm::toString(goal,""));
  
  int numYears = s.fD.finalT;
  double tol = 0.01;

  std::vector<double> scaleD;
  njm::fromFile(scaleD, njm::sett.srcExt("gDistRaw.txt"));
  double pastScale = 1.0;
  double currScale = getDPow(s.modelGen_r.getPar({"power"})[0],
  			     s.modelGen_r.getPar({"alpha"})[0],
  			     s.fD.caves);

  std::vector<double> par = s.modelGen_r.getPar();
  rescaleD(pastScale,currScale,scaleD);
  s.fD.gDist = scaleD;
  s.preCompData();
  s.modelGen_r = MG(s.fD);
  s.modelGen_r.putPar(par.begin());
  s.modelEst_r = MG(s.fD);
  s.modelEst_r.putPar(par.begin());
  
  double val = rn.run(s,nt,numReps,numYears,starts).smean();
  double scale = 1.025, shrink = .9;
  int above = int(val > goal);
  int iter = 0;

  printf("Iter: %05d  >>>  Current value: %08.6f\r",
	 ++iter, val);

  while(std::abs(val - goal) > tol){
    if(val > goal){
      if(!above)
	scale*=shrink;
      
      // std::for_each(par.begin(),par.end(),
      // 		    [&scale](double & x){x*= 1.0 + scale;});
      s.modelGen_r.linScale(1.0 + scale);
      
      above = 1;
    }
    else{
      if(above)
	scale*=shrink;

      // std::for_each(par.begin(),par.end(),
      // 		    [&scale](double & x){x*= 1.0/(1.0 + scale);});
      s.modelGen_r.linScale(1.0/(1.0 + scale));
      
      above = 0;
    }

    par = s.modelGen_r.getPar();
    s.modelEst_r.putPar(par.begin());

    s.revert();

    
    pastScale = currScale;
    currScale = getDPow(s.modelGen_r.getPar({"power"})[0],
			s.modelGen_r.getPar({"alpha"})[0],
			s.fD.caves);

    par = s.modelGen_r.getPar();
    rescaleD(pastScale,currScale,scaleD);
    
    s.fD.gDist = scaleD;
    s.preCompData();
    s.modelGen_r = MG(s.fD);
    s.modelGen_r.putPar(par.begin());
    s.modelEst_r = MG(s.fD);
    s.modelEst_r.putPar(par.begin());


    val = rn.run(s,nt,numReps,numYears,starts).smean();
    printf("Iter: %05d  >>>  Current value: %08.6f\r", ++iter, val);
    fflush(stdout);
  }

  njm::message("Est. goal: " + njm::toString(val,""));

  // write new distance to file
  njm::toFile(njm::toString(s.fD.gDist,"\n",""), njm::sett.srcExt("gDist.txt"),
	      std::ios_base::out);

  njm::message("par: " + njm::toString(s.modelGen_r.getPar()," ",""));

  return(val);
}


template <class S, class MA, class RM, class NT, class RN>
double TuneGenMA(S & s, const int numReps, const Starts & starts){
  NT nt;
  RN rn;
  
  MA ma;
  RM rm;

  int numYears = s.fD.finalT;
  
  double atTrtStart = rn.run(s,nt,numReps,s.fD.trtStart,starts).smean();
  double atFinalT = rn.run(s,nt,numReps,numYears,starts).smean();
  
  double goal = atTrtStart + 0.5*(atFinalT - atTrtStart);
  njm::message("Goal: " + njm::toString(goal,""));
  double tol = 0.01;

  std::vector<double> par;
  double trt = 1.0;
  
  s.modelGen_r.setPar(std::vector<std::string>({"trtAct","trtPre"}),trt);
  par = s.modelGen_r.getPar();
  s.modelEst_r.putPar(par.begin());
  s.revert();
  
  double val = rm.run(s,ma,numReps,numYears,starts).smean();
  double scale = 1.1, shrink = .9;
  int above = int(val > goal);
  int iter = 0;

  printf("Iter: %05d  >>>  Current value: %08.6f\r",
	 ++iter, val);

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
    par = s.modelGen_r.getPar();
    s.modelEst_r.putPar(par.begin());
    s.revert();

    // std::cout << "par: " << njm::toString(par," ","\n");
    // par = s.modelGen.getPar({"trtAct","trtPre"});
    // std::cout << "par: " << njm::toString(par," ","\n");


    val = rm.run(s,ma,numReps,numYears,starts).smean();
    printf("Iter: %05d  >>>  Current value: %08.6f\r", ++iter, val);
    fflush(stdout);
  }

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
  njm::sett.set(argc,argv);

  {
    typedef ModelTimeExpCavesGDist MG;
    typedef MG ME;

    typedef System<MG,ME> S;
    typedef NoTrt<ME> NT;
    typedef ProximalGDistAgent<ME> PA;
    typedef MyopicAgent<ME> MA;

    typedef ToyFeatures5<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef VanillaRunnerNS<S,NT> RN;
    typedef VanillaRunnerNS<S,PA> RP;
    typedef VanillaRunnerNS<S,MA> RM;
    typedef VanillaRunnerNS<S,RA> RR;

    S s;
    s.modelEst_r = s.modelGen_r;
    s.revert();

    njm::resetSeed();
    int numReps = 500;
    Starts starts(numReps,s.fD.numNodes);
    
    PA pa;
    RP rp;

    RA ra;
    RR rr;
    ra.reset();

    njm::message("Tuning Intercept");

    double valNT = TuneGenNT<S,NT,RN,MG>(s,numReps,starts);

    njm::message("Tuning Treatment");

    double valMA = TuneGenMA<S,MA,RM,NT,RN>(s,numReps,starts);

    double valPA = rp.run(s,pa,numReps,s.fD.finalT,starts).smean();

    double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).smean();

    njm::message(" valNT: " + njm::toString(valNT,"") +
		 "\n" +
		 " valPA: " + njm::toString(valPA,"") +
		 "\n" +
		 " valMA: " + njm::toString(valMA,"") +
		 "\n" +
		 " valRA: " + njm::toString(valRA,""));


    s.modelGen_r.save();

    std::vector<double> par = s.modelGen_r.getPar();

    double priorMeanTrt = (s.modelGen_r.getPar({"trtAct"})[0]
			   + s.modelGen_r.getPar({"trtPre"})[0])/2.0;
    priorMeanTrt *= 4.0;

    // // write new distance matrix to file
    // njm::toFile(s.fD.gDist,njm::sett.srcExt("gDist.txt"),
    // 		std::ios_base::out,"\n","");
    // write prior mean of treatment effect
    njm::toFile(priorMeanTrt,njm::sett.srcExt("priorTrtMean.txt"),
		std::ios_base::out);
  }


  
  njm::sett.clean();
  
  return 0;
}
