#include "m1OssOptim.hpp"


M1OsspOptimTunePar::M1OsspOptimTunePar(){
  int N = 1000;
}

std::vector<double> M1OsspOptimTunePar::getPar() const{
  std::vector<double> par(0);
  return par;
}

void M1OsspOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				GravityTimeInfExpCavesModel>,
			 RankAgent<ToyFeatures2<GravityTimeInfExpCavesModel>,
				   GravityTimeInfExpCavesModel>,
			 GravityTimeInfExpCavesModel>;


template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				GravityTimeInfExpModel>,
			 RankAgent<ToyFeatures2<GravityTimeInfExpModel>,
				   GravityTimeInfExpModel>,
			 GravityTimeInfExpModel>;


template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				GravityTimeInfModel>,
			 RankAgent<ToyFeatures2<GravityTimeInfModel>,
				   GravityTimeInfModel>,
			 GravityTimeInfModel>;


template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				GravityModel>,
			 RankAgent<ToyFeatures2<GravityModel>,
				   GravityModel>,
			 GravityModel>;


template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				RangeModel>,
			 RankAgent<ToyFeatures2<RangeModel>,
				   RangeModel>,
			 RangeModel>;

template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				RadiusModel>,
			 RankAgent<ToyFeatures2<RadiusModel>,
				   RadiusModel>,
			 RadiusModel>;

template class M1OsspOptim<System<GravityTimeInfExpCavesModel,
				CaveModel>,
			 RankAgent<ToyFeatures2<CaveModel>,
				   CaveModel>,
			 CaveModel>;


template <class S, class A, class M>
M1OsspOptim<S,A,M>::M1OsspOptim(){
  name = "M1Ossp";
}


template <class S, class A, class M>
void M1OsspOptim<S,A,M>::reset(){
  tp.A = 30;
  tp.B = 1;
}


template <class S, class A, class M>
void M1OsspOptim<S,A,M>
::optim(const S & system,
	A & agent){

  System<M,M> s(system.sD,system.tD,system.fD,system.dD,
		system.modelEst,system.modelEst);

  M1SpOptim<System<M,M>,RankAgent<ToyFeatures2<M>,M>,M> spo;
  RankAgent<ToyFeatures2<M>,M> ra;
  
  spo.optim(s,ra);

  PlainRunner<System<M,M>,A> runner;

  int i,j;
  std::set<std::pair<std::vector<int>,std::vector<int> > > apSet;
  for(i = 0, j = 0; i < 10000 && j < 1000; ++i){
    std::fill(s.tD.a.begin(),s.tD,a.end(),0);
    std::fill(s.tD.p.begin(),s.tD,p.end(),0);
    
    ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
    apSet.insert(std::pair<std::vector<int>, std::vector<int> >
		 (s.tD.a,s.tD.p));
    j = apSet.size();
  }


  

}



template <class S, class A, class M>
void M1OsspOptim<S,A,M>
::tune(const System<M,M> & system,
       A agent){

  std::cout << "thread "
	    << omp_get_thread_num()
	    << " is tuning!!!!!!!" << std::endl;

}
