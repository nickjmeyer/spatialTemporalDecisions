#include "tuneGen2.hpp"


template <class M>
void scaleD(TuneData<M> * const tuneData, const double & pow){
  int i, I = tuneData->s.fD.numNodes * tuneData->s.fD.numNodes;
  for(i = 0; i < I; ++i)
    tuneData->s.fD.gDist.at(i) = std::pow(tuneData->d.at(i),pow);
  
  tuneData->s.preCompData();
  
  tuneData->s.modelGen_r = M(tuneData->s.fD);
  tuneData->s.modelGen_r.setType(INVALID);
  tuneData->s.modelGen_r.putPar(tuneData->par.begin());
}


template <class M>
void tuneSpread(System<M,M> s, const Starts & starts,
		const int numReps,
		double & pow, double & scale,
		const double goal0, const double goal1){
  size_t iter = 0;
  int status;
  gsl_vector *x, *ss;
  int dim=2;

  TuneData<M> tuneData(s,s.modelGen_r.getPar(),s.fD.gDist,
		       starts,numReps,goal0,goal1);

  x = gsl_vector_alloc(dim);
  gsl_vector_set_all(x,0.0);
  ss = gsl_vector_alloc(dim);
  gsl_vector_set_all(ss,0.1);

  gsl_multimin_function minex_func;
  minex_func.n = dim;
  minex_func.f = &tuneSpreadHelper<M>;
  minex_func.params = &tuneData;

  const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer * m = NULL;
  m = gsl_multimin_fminimizer_alloc(T,dim);
  gsl_multimin_fminimizer_set(m,&minex_func,x,ss);

  double curSize;
  double size = 0.001;

  do{
    ++iter;
    status = gsl_multimin_fminimizer_iterate(m);
    if(status)
      break;
    curSize = gsl_multimin_fminimizer_size(m);
    status = gsl_multimin_test_size(curSize,size);
  } while(status == GSL_CONTINUE && iter < 1000);

  pow = std::exp(gsl_vector_get(m->x,0));
  scale = std::exp(gsl_vector_get(m->x,1));


  gsl_multimin_fminimizer_free(m);
  gsl_vector_free(x);
  gsl_vector_free(ss);
}


template <class M>
double tuneSpreadHelper(const gsl_vector * x, void * params){
  typedef System<M,M> S;
  typedef NoTrt<M> NT;
  typedef VanillaRunnerNS<S,NT> RN;
  TuneData<M> * tuneData = static_cast<TuneData<M> *>(params);
  NT nt;
  RN rn;
  
  double curPow = std::exp(gsl_vector_get(x,0));
  double curScale = std::exp(gsl_vector_get(x,1));

  scaleD<M>(tuneData,curPow);
  
  tuneData->s.modelGen.linScale(curScale);

  tuneData->s.modelEst_r = tuneData->s.modelGen_r;

  double val0 = rn.run(tuneData->s,nt,tuneData->numReps,tuneData->s.fD.trtStart,
		       tuneData->starts).smean();

  double val1 = rn.run(tuneData->s,nt,tuneData->numReps,tuneData->s.fD.finalT,
		       tuneData->starts).smean();

  printf("{%5.4f, %5.4f}\r",val0,val1);
  fflush(stdout);

  double ret = std::pow(val0 - tuneData->goal0,2.0) +
    std::pow(val1 - tuneData->goal1,2.0);

  return ret;
}



// void tuneTrt(System & s){
// }



// template <class M>
// double tuneTrtHelper(const gsl_vector * x, void * params){
//   typedef System<M,M> S;
//   typedef Myopic<M> MA;
//   typedef VanillaRunnerNS<S,MA> RM;
//   TuneData<M> * tuneData = static_cast<TuneData<M>*>(params);
//   MA ma;
//   RM rm;
  
//   double curTrt = gsl_vector_get(x,0);
  
//   tuneData->s.modelGen.setPar({"trtPre","trtAct"},curTrt);
//   tuneData->s.modelEst.setPar({"trtPre","trtAct"},curTrt);

//   return rm.run(tuneData->s,ma,tuneData->numReps,tuneData->s.fD.finalT,
// 		tuneData->starts).smean();
// }




int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef ModelTimeExpCavesGDistTrendPowCon MG;
    
    typedef System<MG,MG> S;
    // typedef NoTrt<ME> NT;
    // typedef ProximalGDistAgent<ME> PA;
    // typedef MyopicAgent<ME> MA;

    // typedef ToyFeatures5<ME> F;
    // typedef RankAgent<F,ME> RA;

    // typedef VanillaRunnerNS<S,NT> RN;
    // typedef VanillaRunnerNS<S,PA> RP;
    // typedef VanillaRunnerNS<S,MA> RM;
    // typedef VanillaRunnerNS<S,RA> RR;

    S s;
    s.modelEst_r = s.modelGen_r;
    s.revert();

    int numReps = 500;
    Starts starts(numReps,s.fD.numNodes);

    double pow,scale,goal0,goal1;
    goal0 = 0.3;
    goal1 = 0.7;
    tuneSpread<MG>(s,starts,numReps,pow,scale,goal0,goal1);
    // PA pa;
    // RP rp;

    // RA ra;
    // RR rr;
    // ra.reset();

  //   njm::message("Tuning Intercept");

  //   double valNT = TuneGenNT<S,NT,RN,MG>(s,numReps,starts);

  //   njm::message("Tuning Treatment");

  //   double valMA = TuneGenMA<S,MA,RM,NT,RN>(s,numReps,starts);

  //   double valPA = rp.run(s,pa,numReps,s.fD.finalT,starts).smean();

  //   double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).smean();

  //   njm::message(" valNT: " + njm::toString(valNT,"") +
  // 		 "\n" +
  // 		 " valPA: " + njm::toString(valPA,"") +
  // 		 "\n" +
  // 		 " valMA: " + njm::toString(valMA,"") +
  // 		 "\n" +
  // 		 " valRA: " + njm::toString(valRA,""));


  //   std::vector<double> par = s.modelGen_r.getPar();

  //   double priorMeanTrt = (s.modelGen_r.getPar({"trtAct"})[0]
  // 			   + s.modelGen_r.getPar({"trtPre"})[0])/2.0;
  //   priorMeanTrt *= 4.0;

  //   // // write new distance matrix to file
  //   // njm::toFile(s.fD.gDist,njm::sett.srcExt("gDist.txt"),
  //   // 		std::ios_base::out,"\n","");
  //   // write prior mean of treatment effect
  //   njm::toFile(priorMeanTrt,njm::sett.srcExt("priorTrtMean.txt"),
  // 		std::ios_base::out);
  }
  
  njm::sett.clean();
  
  return 0;
}
