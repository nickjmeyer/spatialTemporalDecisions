#include "m1OsspOptim.hpp"


M1OsspOptimTunePar::M1OsspOptimTunePar(){
  N = 1000;
  B = 100;
  mcReps = 10;
}


std::vector<double> M1OsspOptimTunePar::getPar() const{
  std::vector<double> par(0);
  return par;
}


void M1OsspOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				GravityTimeInfExpCavesModel>,
			   OsspAgent<GravityTimeInfExpCavesModel>,
			   ToyFeatures2<GravityTimeInfExpCavesModel>,
			   GravityTimeInfExpCavesModel>;


template <class S, class A, class F, class M>
M1OsspOptim<S,A,F,M>::M1OsspOptim(){
  name = "M1Ossp";
}


template <class S, class A, class F, class M>
void M1OsspOptim<S,A,F,M>::reset(){
}


template <class S, class A, class F, class M>
void M1OsspOptim<S,A,F,M>
::optim(const S & system,
	A & agent){

  System<M,M> s(system.sD,system.tD,system.fD,system.dD,
		system.modelEst,system.modelEst);

  M1SpOptim<System<M,M>,RankAgent<F,M>,M> spo;
  spo.tp.tune = 0;
  RankAgent<F,M> ra;
  
  spo.optim(s,ra);


  int i,j,I = 10*tp.N;
  std::set<std::pair<std::vector<int>,std::vector<int> > > apSet;
  for(i = 0, j = 0; i < I && j < tp.N; ++i){
    // zero out trt
    std::fill(s.tD.a.begin(),s.tD.a.end(),0);
    std::fill(s.tD.p.begin(),s.tD.p.end(),0);

    // get trt
    ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
    apSet.insert(std::pair<std::vector<int>, std::vector<int> >
		 (s.tD.a,s.tD.p)); // a is first, p is second
    j = apSet.size();
  }
  // zero out trt for good measure
  std::fill(s.tD.a.begin(),s.tD.a.end(),0);
  std::fill(s.tD.p.begin(),s.tD.p.end(),0);


  // clear all current qvalues and trt
  agent.qvalues.clear();
  agent.aCand.clear();
  agent.pCand.clear();

  // setup iterators
  int J = apSet.size();
  std::set<std::pair<std::vector<int>,std::vector<int> > >::iterator it,end;
  end = apSet.end();

  std::cout << "value: " << s.value() << std::endl;
  std::cout << "apSet: " << J << std::endl;


  // initialize runner and original state
  PlainRunner<System<M,M>,RankAgent<F,M> > runner;
  System<M,M> s0(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.modelEst);
  int b;
  double qvalue;
  for(j = 0, it = apSet.begin(); j < J; ++j, ++it){
    printf("\r% 4d",j);
    fflush(stdout);
    qvalue = 0;
    for(b = 0; b < tp.B; ++b){
      // restore original state
      s.sD_r = s0.sD;
      s.dD_r = s0.dD;
      s.modelGen_r = s0.modelGen_r;
      s.modelEst_r = s0.modelEst_r;

      // write new treatments
      s.tD_r.a = (*it).first;
      s.tD_r.p = (*it).second;

      // the above are the revert containers
      s.revert();
      
      // generate next step
      s.nextPoint();

      s.checkPoint();
      
      // estimate V
      // this also is our estimate of Q
      // if we include past rewards it doesn't effect
      // the ranking of current rewards
      qvalue += runner.run(s,ra,tp.mcReps,s.fD.finalT);
    }
    qvalue /= double(tp.B);

    // place qvalue in agent containers
    agent.qvalues.push_back(qvalue);
    agent.aCand.push_back((*it).first);
    agent.pCand.push_back((*it).second);
  }
  printf("\n");
}



template <class S, class A, class F, class M>
void M1OsspOptim<S,A,F,M>
::tune(const System<M,M> & system,
       A agent){

  std::cout << "thread "
	    << omp_get_thread_num()
	    << " is tuning!!!!!!!" << std::endl;

  throw(1);
}
