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

  a=30;
  b=1;

  tune=1;
}

std::vector<double> M1SgdOptimTunePar::getPar() const{
  return std::vector<double> (0);
}

void M1SgdOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1SgdOptim<System<GravityModel,GravityParam,
				 GravityModel,GravityParam>,
			  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
				       GravityModel,GravityParam>,
			  GravityModel,GravityParam>;

template class M1SgdOptim<System<RangeModel,RangeParam,
				 RangeModel,RangeParam>,
			  RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				       RangeModel,RangeParam>,
			  RangeModel,RangeParam>;

template class M1SgdOptim<System<GravityModel,GravityParam,
				 RangeModel,RangeParam>,
			  RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				       RangeModel,RangeParam>,
			  RangeModel,RangeParam>;

template class M1SgdOptim<System<GravityModel,GravityParam,
				 CaveModel,CaveParam>,
			  RankToyAgent<ToyFeatures2<CaveModel,CaveParam>,
				       CaveModel,CaveParam>,
			  CaveModel,CaveParam>;

template class M1SgdOptim<System<CaveModel,CaveParam,
				 CaveModel,CaveParam>,
			  RankToyAgent<ToyFeatures2<CaveModel,CaveParam>,
				       CaveModel,CaveParam>,
			  CaveModel,CaveParam>;



template <class S, class A, class M, class MP>
M1SgdOptim<S,A,M,MP>::M1SgdOptim(){
  name = "M1Sgd";
}

template <class S, class A, class M, class MP>
void M1SgdOptim<S,A,M,MP>
::optim(const S & system,
	A & agent){

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);

  if(tp.tune == 1 && system.sD.time == (system.fD.trtStart + system.fD.period))
    tune(s,agent);

  printf("[optimize]\n");
  
  PlainRunner<System<M,MP,M,MP>,A> runner;

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
  prev = runner.runEx(s,agent,tp.mcReps,s.fD.finalT);
  
  tp.rate = tp.a/tp.b;

  int iter=0;
  while(!converged){

    for(i=0; i<numPar; i++){
      parAdd.at(i) = tp.jitter*njm::rnorm01();
      parJit.at(i) = par.at(i) + parAdd.at(i);
    }

    
    agent.tp.putPar(par);
    prev = runner.runEx(s,agent,tp.mcReps,s.fD.finalT);    

    agent.tp.putPar(parJit);
    curr = runner.runEx(s,agent,tp.mcReps,s.fD.finalT);    

    // if(omp_get_thread_num() == 0)
    //   std::cout << "iter: " + njm::toString(iter,"",4,0) +
    // 	" || " + njm::toString(prev.first,"",6,4) + " - " +
    // 	njm::toString(curr.first,"",6,4) + " -> " +
    // 	njm::toString(curr.second,"",6,4) + " || " +
    // 	njm::toString(par,", ","") << "\r" << std::flush;
    

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



template <class S, class A, class M, class MP>
void M1SgdOptim<S,A,M,MP>
::tune(const System<M,MP,M,MP> & system,
       A agent){

  printf("[tuning]\n");
  // std::cout << "thread "
  // 	    << omp_get_thread_num()
  // 	    << " is tuning!!!!!!!" << std::endl;
  System<M,MP,M,MP> s(system.sD_r,system.tD_r,system.fD,system.dD_r,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);
  s.modelEst.fitType = MLE;
  s.fD.finalT = s.sD.time + 2*s.fD.period;

  M1SgdOptim<System<M,MP,M,MP>,A,M,MP> o;
  o.tp.tune = 0;
  
  TuneRunner<System<M,MP,M,MP>,A,
	     M1SgdOptim<System<M,MP,M,MP>,A,M,MP> > r;
  
  std::vector<double> scale;
  scale.push_back(0.5);
  scale.push_back(1.0);
  scale.push_back(2.0);

  
  int i,j;
  std::vector<std::pair<double,double> > abVals;
  for(i = 0; i < (int)scale.size(); ++i)
    for(j = 0; j < (int)scale.size(); ++j)
      abVals.push_back(std::pair<double,double>(30*scale.at(i),
						1*scale.at(j)));

  int numAbVals=abVals.size();
  double val,minVal=1.0,bestA=10,bestB=100;
  for(i=0; i<numAbVals; i++){
    o.tp.a=abVals.at(i).first;
    o.tp.b=abVals.at(i).second;

    printf("[setting](% 4d)\n",i);
    val = r.run(s,agent,o,1,s.fD.finalT);
    if(val < minVal){
      bestA = o.tp.a;
      bestB = o.tp.b;

      minVal = val;
    }
  }

  tp.a = bestA;
  tp.b = bestB;
}
