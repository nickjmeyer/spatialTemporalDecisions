#include "m1HybridOptim.hpp"


M1HybridOptimTunePar::M1HybridOptimTunePar(){
  topWeights = 1;
  
  // for simple optim
  mcRepsSimple = 250;

  // for sgd optim
  mcRepsSgd = 100;
  aSgd = 30;
  bSgd = 1;
  tolSgd = .0001;
}

std::vector<double> M1HybridOptimTunePar::getPar() const{
  return std::vector<double> (0);
}

void M1HybridOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1HybridOptim<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     RankToyAgent<ToyFeatures0<GravityModel,
						       GravityParam>,
					  GravityModel,GravityParam>,
			     GravityModel,GravityParam>;
template class M1HybridOptim<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     RankToyAgent<ToyFeatures1<GravityModel,
						       GravityParam>,
					  GravityModel,GravityParam>,
			     GravityModel,GravityParam>;

template class M1HybridOptim<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     RankToyAgent<ToyFeatures2<GravityModel,
						       GravityParam>,
					  GravityModel,GravityParam>,
			     GravityModel,GravityParam>;

template class M1HybridOptim<System<RangeModel,RangeParam,
				    RangeModel,RangeParam>,
			     RankToyAgent<ToyFeatures2<RangeModel,
						       RangeParam>,
					  RangeModel,RangeParam>,
			     RangeModel,RangeParam>;

template class M1HybridOptim<System<GravityModel,GravityParam,
				    RangeModel,RangeParam>,
			     RankToyAgent<ToyFeatures2<RangeModel,
						       RangeParam>,
					  RangeModel,RangeParam>,
			     RangeModel,RangeParam>;


template class M1HybridOptim<System<EbolaModel,EbolaParam,
				    EbolaModel,EbolaParam>,
			     RankToyAgent<ToyFeatures1<EbolaModel,
						       EbolaParam>,
					  EbolaModel,EbolaParam>,
			     EbolaModel,EbolaParam>;


template <class S, class A, class M, class MP>
M1HybridOptim<S,A,M,MP>::M1HybridOptim(){
  name = "M1Hybrid";
}

template <class S, class A, class M, class MP>
void M1HybridOptim<S,A,M,MP>
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


      // double j
      w.at(j)=2;
      w.at(k)=1;
      weights.push_back(w);

      // double k
      w.at(j)=1;
      w.at(k)=2;
      weights.push_back(w);

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
    p.push(std::pair<double,int>(-runner.run(s,agent,tp.mcRepsSimple,
					     s.fD.finalT),
				 i));
  }

  // create the optimzation object and set tuning parameters
  M1SgdOptim<System<M,MP,M,MP>,A,M,MP> sO;
  sO.tp.a = tp.aSgd;
  sO.tp.b = tp.bSgd;
  sO.tp.mcReps = tp.mcRepsSgd;
  sO.tp.tol = tp.tolSgd;

  // optimize the top sets of weights
  std::priority_queue< std::pair<double,int> > pC=p; // pC is a copy
  int NO = N + std::min(N,tp.topWeights);
  int top;
  for(i=N; i<NO; i++){
    top = pC.top().second; // which index is next
    pC.pop(); // pop the top
    
    agent.tp.putPar(weights.at(top)); // assign weights
    sO.optim(s,agent); // optimize weights
    
    weights.push_back(agent.tp.getPar()); // add weights to vector
    
    p.push(std::pair<double,int>(-runner.run(s,agent,tp.mcRepsSimple,
					     s.fD.finalT),
				 i)); // i starts at N
  }

  // assign optimized par to the agent
  agent.tp.putPar(weights.at(p.top().second));
}


