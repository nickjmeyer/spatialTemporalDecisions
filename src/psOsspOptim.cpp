#include "psOsspOptim.hpp"


PsOsspOptimTunePar::PsOsspOptimTunePar(){
  // N = 1000;
  N = 100;
  B = 100;
  mcReps = 10;

  corrGoal = 0.8;
  
  // N = 10;
  // B = 5;
  // mcReps = 2;
}


std::vector<double> PsOsspOptimTunePar::getPar() const{
  std::vector<double> par(0);
  return par;
}


void PsOsspOptimTunePar::putPar(const std::vector<double> & par){
}


template class
PsOsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
	    ModelTimeExpCavesGDistTrendPowCon>;



template <class S, class A, class M>
PsOsspOptim<S,A,M>::PsOsspOptim(){
  name = "PsOssp";
}


template <class S, class A, class M>
void PsOsspOptim<S,A,M>::reset(){
}


template <class S, class A, class M>
void PsOsspOptim<S,A,M>
::optim(const S & system,
	A & agent){

  System<M,M> s(system.sD,system.tD,system.fD,system.dD,
		system.modelEst,system.modelEst);

  ProxStocGDistAgent<M> ps;
  ps.tp.corrGoal = tp.corrGoal;


  int i,j,I = 10*tp.N;
  std::set<std::pair<std::vector<int>,std::vector<int> > > apSet;
  for(i = 0, j = 0; i < I && j < tp.N; ++i){
    // zero out trt
    std::fill(s.tD.a.begin(),s.tD.a.end(),0);
    std::fill(s.tD.p.begin(),s.tD.p.end(),0);

    // get trt
    ps.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
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


  // initialize runner and original state
  PlainRunner<System<M,M>, ProximalGDistAgent<M> > runner;
  ProximalGDistAgent<M> pa;
  System<M,M> s0(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.modelEst);
  int b;
  double qvalue;
  for(j = 0, it = apSet.begin(); j < J; ++j, ++it){
    // printf("\r% 4d",j);
    // fflush(stdout);
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
      qvalue += - runner.run(s,pa,tp.mcReps,s.fD.finalT).smean();
    }
    qvalue /= double(tp.B);

    // place qvalue in agent containers
    agent.qvalues.push_back(qvalue);
    agent.aCand.push_back((*it).first);
    agent.pCand.push_back((*it).second);
  }
  // printf("\n");
}



