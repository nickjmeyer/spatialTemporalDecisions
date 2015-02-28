#include "m1SimpleOptim.hpp"


M1SimpleOptimTunePar::M1SimpleOptimTunePar(){
  mcReps = 250;
}

std::vector<double> M1SimpleOptimTunePar::getPar() const{
  return std::vector<double> (0);
}

void M1SimpleOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1SimpleOptim<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     RankToyAgent<ToyFeatures2<GravityModel,
						       GravityParam>,
					  GravityModel,GravityParam>,
			     GravityModel,GravityParam>;

template class M1SimpleOptim<System<RangeModel,RangeParam,
				    RangeModel,RangeParam>,
			     RankToyAgent<ToyFeatures2<RangeModel,
						       RangeParam>,
					  RangeModel,RangeParam>,
			     RangeModel,RangeParam>;

template class M1SimpleOptim<System<GravityModel,GravityParam,
				    RangeModel,RangeParam>,
			     RankToyAgent<ToyFeatures2<RangeModel,
						       RangeParam>,
					  RangeModel,RangeParam>,
			     RangeModel,RangeParam>;

template class M1SimpleOptim<System<GravityModel,GravityParam,
				    CaveModel,CaveParam>,
			     RankToyAgent<ToyFeatures2<CaveModel,
						       CaveParam>,
					  CaveModel,CaveParam>,
			     CaveModel,CaveParam>;


template <class S, class A, class M , class MP>
M1SimpleOptim<S,A,M,MP>::M1SimpleOptim(){
  name = "M1Simple";
}

template <class S, class A, class M, class MP>
void M1SimpleOptim<S,A,M,MP>
::optim(const S & system,
	A & agent){

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);
  
  PlainRunner<System<M,MP,M,MP>,A> runner;

  int i,j,k;
  std::vector< std::vector<double> > weights;
  std::vector<double> w;
  w.resize(agent.f.numFeatures);
  std::fill(w.begin(),w.end(),0);
  for(j=0; j<agent.f.numFeatures; j++){
    std::fill(w.begin(),w.end(),0);
    w.at(j)=1;
    weights.push_back(w);
    
    w.at(j)=-1;
    weights.push_back(w);
    
    for(k=(j+1); k<agent.f.numFeatures; k++){
      std::fill(w.begin(),w.end(),0);
      w.at(j)=1;
      w.at(k)=1;
      weights.push_back(w);

      // w.at(j)=-1;
      // w.at(k)=1;
      // weights.push_back(w);

      // w.at(j)=1;
      // w.at(k)=-1;
      // weights.push_back(w);

      // w.at(j)=-1;
      // w.at(k)=-1;
      // weights.push_back(w);

      // double j
      w.at(j)=2;
      w.at(k)=1;
      weights.push_back(w);

      // w.at(j)=-2;
      // w.at(k)=1;
      // weights.push_back(w);

      // w.at(j)=2;
      // w.at(k)=-1;
      // weights.push_back(w);

      // w.at(j)=-2;
      // w.at(k)=-1;
      // weights.push_back(w);

      // double k
      w.at(j)=1;
      w.at(k)=2;
      weights.push_back(w);

      // w.at(j)=-1;
      // w.at(k)=2;
      // weights.push_back(w);

      // w.at(j)=1;
      // w.at(k)=-2;
      // weights.push_back(w);

      // w.at(j)=-1;
      // w.at(k)=-2;
      // weights.push_back(w);
    }
  }
  std::fill(w.begin(),w.end(),1);
  weights.push_back(w);
  
  std::fill(w.begin(),w.end(),-1);
  weights.push_back(w);
  

  int N=weights.size();
  std::priority_queue< std::pair<double,int> > p;
  for(i=0; i<N; i++){
    agent.tp.putPar(weights.at(i));
    p.push(std::pair<double,int>(-runner.run(s,agent,tp.mcReps,
					     s.fD.finalT),
				 i));
  }

  // assign optimized par to the agent
  agent.tp.putPar(weights.at(p.top().second));
}


