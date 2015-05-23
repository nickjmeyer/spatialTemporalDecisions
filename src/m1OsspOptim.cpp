#include "m1OsspOptim.hpp"


M1OsspOptimTunePar::M1OsspOptimTunePar(){
  // N = 1000;
  N = 100;
  B = 100;
  mcReps = 10;

  jitterScale = 4.0;
  
  // N = 10;
  // B = 5;
  // mcReps = 2;
}


std::vector<double> M1OsspOptimTunePar::getPar() const{
  std::vector<double> par(0);
  return par;
}


void M1OsspOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   OsspAgent<ModelTimeExpCavesGDist>,
			   ToyFeatures0<ModelTimeExpCavesGDist>,
			   ModelTimeExpCavesGDist>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelGDist>,
			   OsspAgent<ModelGDist>,
			   ToyFeatures0<ModelGDist>,
			   ModelGDist>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   OsspAgent<ModelTimeExpCavesGDist>,
			   ToyFeatures4<ModelTimeExpCavesGDist>,
			   ModelTimeExpCavesGDist>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelRadius>,
			   OsspAgent<ModelRadius>,
			   ToyFeatures4<ModelRadius>,
			   ModelRadius>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelGDist>,
			   OsspAgent<ModelGDist>,
			   ToyFeatures4<ModelGDist>,
			   ModelGDist>;


template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   OsspAgent<ModelTimeExpCavesGDist>,
			   ToyFeatures5<ModelTimeExpCavesGDist>,
			   ModelTimeExpCavesGDist>;

template class M1OsspOptim<System<ModelTimeGDistTrendPow,
				  ModelTimeGDistTrendPow>,
			   OsspAgent<ModelTimeGDistTrendPow>,
			   ToyFeatures5<ModelTimeGDistTrendPow>,
			   ModelTimeGDistTrendPow>;

template class M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
				  ModelTimeExpCavesGDistTrendPowCon>,
			   OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
			   ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
			   ModelTimeExpCavesGDistTrendPowCon>;

template class M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
				  ModelTimeExpCavesEDist>,
			   OsspAgent<ModelTimeExpCavesEDist>,
			   ToyFeatures5<ModelTimeExpCavesEDist>,
			   ModelTimeExpCavesEDist>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesEDist>,
			   OsspAgent<ModelTimeExpCavesEDist>,
			   ToyFeatures5<ModelTimeExpCavesEDist>,
			   ModelTimeExpCavesEDist>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelRadius>,
			   OsspAgent<ModelRadius>,
			   ToyFeatures5<ModelRadius>,
			   ModelRadius>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelGDist>,
			   OsspAgent<ModelGDist>,
			   ToyFeatures5<ModelGDist>,
			   ModelGDist>;


template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelCovar>,
			   OsspAgent<ModelCovar>,
			   ToyFeatures5<ModelCovar>,
			   ModelCovar>;



template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   OsspAgent<ModelTimeExpCavesGDist>,
			   WnsFeatures1<ModelTimeExpCavesGDist>,
			   ModelTimeExpCavesGDist>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelRadius>,
			   OsspAgent<ModelRadius>,
			   WnsFeatures1<ModelRadius>,
			   ModelRadius>;

template class M1OsspOptim<System<ModelTimeExpCavesGDist,
				  ModelGDist>,
			   OsspAgent<ModelGDist>,
			   WnsFeatures1<ModelGDist>,
			   ModelGDist>;



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
  RankAgent<F,M> ra;
  ra.tp.jitterScale = tp.jitterScale;
  
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


  // initialize runner and original state
  PlainRunner<System<M,M>,RankAgent<F,M> > runner;
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
      qvalue += - runner.run(s,ra,tp.mcReps,s.fD.finalT).smean();
    }
    qvalue /= double(tp.B);

    // place qvalue in agent containers
    agent.qvalues.push_back(qvalue);
    agent.aCand.push_back((*it).first);
    agent.pCand.push_back((*it).second);
  }
  // printf("\n");
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
