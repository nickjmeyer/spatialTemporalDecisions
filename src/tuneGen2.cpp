#include "tuneGen2.hpp"

template<class M>
class TuneData {
 public:
  System<M,M> s;
  Starts starts;
  int numReps;
};

void tuneSpread(System & s);

void tuneSpread(System & s, const Starts & starts,
		const int numReps){
  size_t iter = 0;
  int status;
  gsl_vector *x, *ss;
  int i,dim=2;

  TuneData tuneData(s,starts,numReps);

  x = gsl_vector_alloc(dim);
  gsl_vector_set_all(x,1.0);
  ss = gsl_vector_alloc(dim);
  gsl_vector_set_all(ss,0.1);

  gsl_multimin_function minex_func;
  minex_func.n = dim;
  minex_func.f = &tuneSpreadHelper;
  minex_func.params = & tuneData;

  const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer * s = NULL;
  s = gsl_multimin_fminmizer_alloc(T,dim);
  gsl_multimin_fminimizer_set(s,&minex_func,x,ss);

  double curSize;
  double size = 0.001;

  do{
    ++iter;
    status = gsl_multimin_fminimizer_iterate(s);
    if(status)
      break;
    curSize = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(curSize,size);
  } while(status == GSL_CONTINUE && iter < 1000);



  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(ss);
}

double tuneSpreadHelper(const gsl_vector * x, void * params);

template <class M>
double tuneSpreadHelper(const gsl_vector * x, void * params){
  typedef System<M,M> S;
  typedef NoTrt<M> NT;
  typedef VanillaRunnerNS<S,NT> RN;
  TuneData<M> * tuneData = static_cast<TuneData<M>*>(params);
  NT nt;
  RN rn;
  
  static double pastScale = 1.0;

  double curPow = std::exp(gsl_vector_get(x,0));
  double curScale = gsl_vector_get(x,1);
  
  tuneData->s.modelGen.setPar("dPow",curPow);
  tuneData->s.modelGen.linScale(gsl_vector_get(x,1)/pastScale);
  tuneData->s.modelEst.setPar("dPow",curPow);
  tuneData->s.modelEst.linScale(gsl_vector_get(x,1)/pastScale);
  pastScale = gsl_vector_get(x,1);

  return rn.run(tuneData->s,nt,tuneData->numReps,tuneData->s.fD.finalT,
		tuneData->starts).smean();
}



void tuneTrt(System & s);

void tuneTrt(System & s){
}

double tuneTrtHelper(const gsl_vector * x, void * params);

template <class M>
double tuneTrtHelper(const gsl_vector * x, void * params){
  typedef System<M,M> S;
  typedef Myopic<M> MA;
  typedef VanillaRunnerNS<S,MA> RM;
  TuneData<M> * tuneData = static_cast<TuneData<M>*>(params);
  MA ma;
  RM rm;
  
  double curTrt = gsl_vector_get(x,0);
  
  tuneData->s.modelGen.setPar({"trtPre","trtAct"},curTrt);
  tuneData->s.modelEst.setPar({"trtPre","trtAct"},curTrt);

  return rm.run(tuneData->s,ma,tuneData->numReps,tuneData->s.fD.finalT,
		tuneData->starts).smean();
}




int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef ModelTimeExpCavesGDistTrendPowCon MG;
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

    int numReps = 500;
    Starts starts(numReps,s.fD.numNodes);
    
    PA pa;
    RP rp;

    RA ra;
    RR rr;
    // ra.reset();

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
