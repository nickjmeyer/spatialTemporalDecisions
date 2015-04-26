#include "tuneGenWNS.hpp"

template <class M>
double getAlpha(const M & m,
		const FixedData & fD){
  std::vector<double> par = m.getPar();
  return par.at(1 + fD.numCovar);
}

template <class M>
void putAlpha(const double & alpha,
	      M & m,
	      const FixedData & fD){
  std::vector<double> par = m.getPar();
  par.at(1 + fD.numCovar) = alpha;
  m.putPar(par.begin());
}

template <class M>
double getPower(const M & m,
		const FixedData & fD){
  std::vector<double> par = m.getPar();
  return par.at(1 + fD.numCovar + 1);
}

template <class M>
void putPower(const double & power,
	      M & m,
	      const FixedData & fD){
  std::vector<double> par = m.getPar();
  par.at(1 + fD.numCovar + 1) = power;
  m.putPar(par.begin());
}

template <class M>
double getActTrt(const M & m,
		 const FixedData & fD){
  std::vector<double> par = m.getPar();
  return par.at(par.size() - 2);
}

template <class M>
void putActTrt(const double & trt,
	       M & m,
	       const FixedData & fD){
  std::vector<double> par = m.getPar();
  par.at(par.size() - 2) = trt;
  m.putPar(par.begin());
}


template <class M>
double getPreTrt(const M & m,
		 const FixedData & fD){
  std::vector<double> par = m.getPar();
  return par.at(par.size() - 1);
}

template <class M>
void putPreTrt(const double & trt,
	       M & m,
	       const FixedData & fD){
  std::vector<double> par = m.getPar();
  par.at(par.size() - 1) = trt;
  m.putPar(par.begin());
}


template <class S, class PA, class RP>
double TuneGenPA(S & s,const Starts & starts){
  double trtSize = s.modelGen.tuneTrt(s.fD);

  putActTrt(trtSize,s.modelGen_r,s.fD);
  putPreTrt(trtSize,s.modelGen_r,s.fD);
  putActTrt(trtSize,s.modelEst_r,s.fD);
  putPreTrt(trtSize,s.modelEst_r,s.fD);

  PA pa;
  RP rp;

  return rp.run(s,pa,500,s.fD.finalT,starts);
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef ModelGravity GM;
    typedef GM EM;

    typedef System<GM,EM> S;
    typedef NoTrt<EM> NT;
    typedef ProximalAgent<EM> PA;
    typedef MyopicAgent<EM> MA;

    typedef WnsFeatures0<EM> F;
    typedef RankAgent<F,EM> RA;

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

    NT nt;
    MA ma;
    RM rm;

    RN rn;
    RA ra;
    RR rr;
    ra.reset();

    njm::message("Tuning Treatment");

    double valNT = rn.run(s,nt,500,s.fD.finalT,starts);

    double valPA = TuneGenPA<S,PA,RP>(s,starts);

    double valMA = rm.run(s,ma,500,s.fD.finalT,starts);

    double valRA = rr.run(s,ra,500,s.fD.finalT,starts);

    njm::message(" valNT: " + njm::toString(valNT,"") +
		 "\n" +
		 " valPA: " + njm::toString(valPA,"") +
		 "\n" +
		 " valMA: " + njm::toString(valMA,"") +
		 "\n" +
		 " valRA: " + njm::toString(valRA,""));

    s.modelGen_r.save();
  }

  
  {
    typedef ModelTimeExpCaves GM;
    typedef GM EM;

    typedef System<GM,EM> S;
    typedef NoTrt<EM> NT;
    typedef ProximalAgent<EM> PA;
    typedef MyopicAgent<EM> MA;

    typedef WnsFeatures0<EM> F;
    typedef RankAgent<F,EM> RA;

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

    NT nt;
    MA ma;
    RM rm;

    RN rn;
    RA ra;
    RR rr;
    ra.reset();

    njm::message("Tuning Treatment");

    double valNT = rn.run(s,nt,500,s.fD.finalT,starts);

    double valPA = TuneGenPA<S,PA,RP>(s,starts);

    double valMA = rm.run(s,ma,500,s.fD.finalT,starts);

    double valRA = rr.run(s,ra,500,s.fD.finalT,starts);

    njm::message(" valNT: " + njm::toString(valNT,"") +
		 "\n" +
		 " valPA: " + njm::toString(valPA,"") +
		 "\n" +
		 " valMA: " + njm::toString(valMA,"") +
		 "\n" +
		 " valRA: " + njm::toString(valRA,""));

    s.modelGen_r.save();

    s.modelGen_r.save();


    double priorMeanTrt = (getPreTrt(s.modelGen_r,s.fD)
			   + getPreTrt(s.modelGen_r,s.fD))/2.0;
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
