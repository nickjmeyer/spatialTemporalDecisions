#include "m1SgdOptim.hpp"


M1SgdOptimTunePar::M1SgdOptimTunePar(){
  jitter = .05;
  // mcReps = 10;
  mcReps = 100;
  // tol = .025;
  
  tol = .0001;
  rate = 20;
  rateDecay = .975;

  momRate=.5;

  a=5;
  b=5;
}

std::vector<double> M1SgdOptimTunePar::getPar() const{
  return std::vector<double> (0);
}

void M1SgdOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1SgdOptim<System<GravityModel,GravityParam,
				 GravityModel,GravityParam>,
			  RankToyAgent<ToyFeatures0<GravityModel,GravityParam>,
				       GravityModel,GravityParam> >;
template class M1SgdOptim<System<GravityModel,GravityParam,
				 GravityModel,GravityParam>,
			  RankToyAgent<ToyFeatures1<GravityModel,GravityParam>,
				       GravityModel,GravityParam> >;
template class M1SgdOptim<System<GravityModel,GravityParam,
				 GravityModel,GravityParam>,
			  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
				       GravityModel,GravityParam> >;

template class M1SgdOptim<System<GravityModel,GravityParam,
				 RangeModel,RangeParam>,
			  RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				       RangeModel,RangeParam> >;

template class M1SgdOptim<System<EbolaModel,EbolaParam,
				 EbolaModel,EbolaParam>,
			  RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
				       EbolaModel,EbolaParam> >;


template <class System, class Agent>
M1SgdOptim<System,Agent>::M1SgdOptim(){
  name = "M1Sgd";
}

template <class System, class Agent>
void M1SgdOptim<System,Agent>
::optim(System system,
	Agent & agent){
  
  PlainRunner<System,Agent> runner;

  system.checkPoint();
  
  std::vector<double> par=agent.tp.getPar();
  int i,sameRep=0,converged=0,numPar = par.size();
  std::vector<double> parJit(numPar);
  std::vector<double> parAdd(numPar);
  std::vector<double> parNew(numPar);
  std::vector<double> mom(numPar);
  std::fill(parJit.begin(),parJit.end(),0.0);
  std::fill(parAdd.begin(),parAdd.end(),0.0);
  std::fill(parNew.begin(),parNew.end(),0.0);
  std::fill(mom.begin(),mom.end(),0.0);


  std::pair<double,double> curr,prev;
  prev = runner.runEx(system,agent,tp.mcReps,system.fD.finalT);
  
  tp.rate = tp.a/tp.b;

  int iter=0;
  while(!converged){

    for(i=0; i<numPar; i++){
      parAdd.at(i) = tp.jitter*njm::rnorm01();
      parJit.at(i) = par.at(i) + parAdd.at(i);
    }

    
    agent.tp.putPar(par);
    prev = runner.runEx(system,agent,tp.mcReps,system.fD.finalT);    

    agent.tp.putPar(parJit);
    curr = runner.runEx(system,agent,tp.mcReps,system.fD.finalT);    

    // njm::message("iter: " + njm::toString(iter,"",4,0) +
    // 		 " || " + njm::toString(prev.first,"",6,4) + " - " +
    // 		 njm::toString(curr.first,"",6,4) + " -> " +
    // 		 njm::toString(curr.second,"",6,4) + " || " +
    // 		 njm::toString(par,", ",""));
    

    for(i=0; i<numPar; i++)
      mom.at(i) = tp.momRate * mom.at(i) +
	tp.rate * (curr.first - prev.first) * parAdd.at(i);
      
    for(i=0; i<numPar; i++)
      parNew.at(i) = par.at(i) - mom.at(i);

    njm::l2norm(parNew); // normalize
    
    if(njm::l2norm(parNew,par) < tp.tol){ // check if difference is small
      sameRep++;
      if(sameRep < 0)
	sameRep = 0;
      else if(sameRep == 3) // if small difference 3 times, it converged
	converged = 1;
    }

    iter++;
    
    tp.rate = tp.a/(tp.b + iter);
    // tp.rate *= tp.rateDecay;
    par = parNew;
  }

  agent.tp.putPar(par); // assign optimized par to the agent
}



template <class System, class Agent>
void M1SgdOptim<System,Agent>
::tune(System system,
       Agent & agent){

  std::vector<double> aVals;
  aVals.push_back(1);
  aVals.push_back(3);
  aVals.push_back(7);
  aVals.push_back(10);
  aVals.push_back(13);
  aVals.push_back(17);
  aVals.push_back(20);
  aVals.push_back(40);
  
  int i;
  std::vector<std::pair<double,double> > abVals;
  for(i=0; i<(int)aVals.size(); i++)
    abVals.push_back(std::pair<double,double>(aVals.at(i),1));
  abVals.push_back(std::pair<double,double>(1,40));
  abVals.push_back(std::pair<double,double>(40,40));

  PlainRunner<System,Agent> pR;
  
  int numAbVals=abVals.size();
  double val,minVal=1.0,bestA=40,bestB=1;
  for(i=0; i<numAbVals; i++){
    tp.a=abVals.at(i).first;
    tp.b=abVals.at(i).second;

    val=pR.run(system,agent,150,system.fD.finalT);
    if(val<=minVal){
      bestA = abVals.at(i).first;
      bestB = abVals.at(i).second;

      minVal = val;
    }
  }

  tp.a = bestA;
  tp.b = bestB;
}
