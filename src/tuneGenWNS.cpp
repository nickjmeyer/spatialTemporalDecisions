#include "tuneGenWNS.hpp"

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
    typedef MyopicAgent<EM,EP> MA;

    typedef ToyFeatures2<EM,EP> F;
    typedef RankAgent<F,EM,EP> RA;

    typedef VanillaRunnerNS<S,NT> RN;
    typedef VanillaRunnerNS<S,PA> RP;
    typedef VanillaRunnerNS<S,MA> RM;
    typedef VanillaRunnerNS<S,RA> RR;

    S s;
    s.paramEst_r = s.paramGen_r;
    s.reset();

    NT nt;
    MA ma;
    RM rm;

    RN rn;
    RA ra;
    RR rr;
    ra.reset();

    njm::message("Tuning Treatment");

    double valNT = rn.run(s,nt,500,s.fD.finalT);

    double valPA = TuneGenPA<S,PA,RP>(s);

    double valMA = rm.run(s,ma,500,s.fD.finalT);

    double valRA = rr.run(s,ra,500,s.fD.finalT);

    njm::message(" intcp: " + njm::toString(s.paramGen_r.intcp,"") +
		 "\n" +
		 " alpha: " + njm::toString(s.paramGen_r.alpha,"") +
		 "\n" +
		 " power: " + njm::toString(s.paramGen_r.power,"") +
		 "\n" +
		 "trtPre: " + njm::toString(s.paramGen_r.trtPre,"") +
		 "\n" +
		 "trtAct: " + njm::toString(s.paramGen_r.trtAct,"") +
		 "\n\n" +
		 " valNT: " + njm::toString(valNT,"") +
		 "\n" +
		 " valPA: " + njm::toString(valPA,"") +
		 "\n" +
		 " valMA: " + njm::toString(valMA,"") +
		 "\n" +
		 " valRA: " + njm::toString(valRA,""));

    s.paramGen_r.save();
  }

  
  {
    typedef GravityTimeInfExpCavesModel GM;
    typedef GravityTimeInfExpCavesParam GP;
    typedef GM EM;
    typedef GP EP;

    typedef System<GM,GP,EM,EP> S;
    typedef NoTrt<EM,EP> NT;
    typedef ProximalAgent<EM,EP> PA;
    typedef MyopicAgent<EM,EP> MA;

    typedef ToyFeatures2<EM,EP> F;
    typedef RankAgent<F,EM,EP> RA;

    typedef VanillaRunnerNS<S,NT> RN;
    typedef VanillaRunnerNS<S,PA> RP;
    typedef VanillaRunnerNS<S,MA> RM;
    typedef VanillaRunnerNS<S,RA> RR;

    S s;
    s.paramEst_r = s.paramGen_r;
    s.reset();

    NT nt;
    MA ma;
    RM rm;

    RN rn;
    RA ra;
    RR rr;
    ra.reset();

    njm::message("Tuning Treatment");

    double valNT = rn.run(s,nt,500,s.fD.finalT);

    double valPA = TuneGenPA<S,PA,RP>(s);

    double valMA = rm.run(s,ma,500,s.fD.finalT);

    double valRA = rr.run(s,ra,500,s.fD.finalT);

    njm::message(" intcp: " + njm::toString(s.paramGen_r.intcp,"") +
		 "\n" +
		 " alpha: " + njm::toString(s.paramGen_r.alpha,"") +
		 "\n" +
		 " power: " + njm::toString(s.paramGen_r.power,"") +
		 "\n" +
		 "trtPre: " + njm::toString(s.paramGen_r.trtPre,"") +
		 "\n" +
		 "trtAct: " + njm::toString(s.paramGen_r.trtAct,"") +
		 "\n\n" +
		 " valNT: " + njm::toString(valNT,"") +
		 "\n" +
		 " valPA: " + njm::toString(valPA,"") +
		 "\n" +
		 " valMA: " + njm::toString(valMA,"") +
		 "\n" +
		 " valRA: " + njm::toString(valRA,""));

    s.paramGen_r.save();

    s.paramGen_r.save();


    double priorMeanTrt = (s.paramGen_r.trtPre + s.paramGen_r.trtAct)/2.0;
    priorMeanTrt *= 4.0;

    // write new distance matrix to file
    njm::toFile(s.fD.dist,njm::sett.srcExt("d.txt"),
		std::ios_base::out,"\n","");
    // write prior mean of treatment effect
    njm::toFile(priorMeanTrt,njm::sett.srcExt("priorTrtMean.txt"),
		std::ios_base::out);
  }


  
  njm::sett.clean();
  
  return 0;
}
