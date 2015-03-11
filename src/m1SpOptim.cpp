#include "m1SpOptim.hpp"


M1SpOptimTunePar::M1SpOptimTunePar(){
  mcReps = 100;

  C = 2.0;

  t = 2.0;

  ell = 1.25;

  muMin = 0.1;

  A = 50;
  B = 10;

  tune = 1;
}

std::vector<double> M1SpOptimTunePar::getPar() const{
  std::vector<double> par = {A,B};
  return par;
}

void M1SpOptimTunePar::putPar(const std::vector<double> & par){
  std::vector<double>::const_iterator it;
  it = par.begin();
  A = *it++;
  B = *it++;
}


template class M1SpOptim<System<GravityModel,GravityParam,
				GravityModel,GravityParam>,
			 RankAgent<ToyFeatures2<GravityModel,GravityParam>,
				   GravityModel,GravityParam>,
			 GravityModel,GravityParam>;

template class M1SpOptim<System<GravityTimeInfModel,GravityTimeInfParam,
				GravityTimeInfModel,GravityTimeInfParam>,
			 RankAgent<ToyFeatures2<GravityTimeInfModel,
						GravityTimeInfParam>,
				   GravityTimeInfModel,GravityTimeInfParam>,
			 GravityTimeInfModel,GravityTimeInfParam>;

template class M1SpOptim<System<GravityModel,GravityParam,
				RangeModel,RangeParam>,
			 RankAgent<ToyFeatures2<RangeModel,RangeParam>,
				   RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class M1SpOptim<System<GravityModel,GravityParam,
				CaveModel,CaveParam>,
			 RankAgent<ToyFeatures2<CaveModel,CaveParam>,
				   CaveModel,CaveParam>,
			 CaveModel,CaveParam>;

template class M1SpOptim<System<RangeModel,RangeParam,
				RangeModel,RangeParam>,
			 RankAgent<ToyFeatures2<RangeModel,RangeParam>,
				   RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class M1SpOptim<System<CaveModel,CaveParam,
				CaveModel,CaveParam>,
			 RankAgent<ToyFeatures2<CaveModel,CaveParam>,
				   CaveModel,CaveParam>,
			 CaveModel,CaveParam>;


template <class S, class A, class M, class MP>
M1SpOptim<S,A,M,MP>::M1SpOptim(){
  name = "M1Sp";
}


template <class S, class A, class M, class MP>
void M1SpOptim<S,A,M,MP>::reset(){
  tp.A = 50;
  tp.B = 10;
}


template <class S, class A, class M, class MP>
void M1SpOptim<S,A,M,MP>
::optim(const S & system,
	A & agent){

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);

  if(tp.tune == 1 && system.sD.time == (system.fD.trtStart + 1))
    tune(s,agent);

  
  PlainRunner<System<M,MP,M,MP>,A> runner;

  std::vector<double> par=agent.tp.getPar();
  int i,converged=0,numPar = par.size();
  std::vector<double> h(numPar,0.0);
  std::vector<double> parPH(numPar,0.0),parMH(numPar,0.0);


  double valP, valM;
  
  int iter=1;
  
  double mu = tp.A/std::pow(tp.B+iter,tp.ell);
  double cm = tp.C / std::pow(iter,tp.t);

  while(!converged){

    for(i=0; i<numPar; i++){
      h.at(i) = (njm::rber(0.5) == 1 ? 1.0 : -1.0) * cm;
      parPH.at(i) = par.at(i) + h.at(i);
      parMH.at(i) = par.at(i) - h.at(i);
    }


    
    agent.tp.putPar(parPH);
    valP = runner.run(s,agent,tp.mcReps,s.fD.finalT);

    agent.tp.putPar(parMH);
    valM = runner.run(s,agent,tp.mcReps,s.fD.finalT);

    
    for(i=0; i<numPar; i++)
      par.at(i) = par.at(i) - mu*(valP - valM)/(2.0*h.at(i));

    
    // if(omp_get_thread_num() == 0)
    //   std::cout << "iter: " + njm::toString(iter,"",4,0) +
    // 	" || " + njm::toString(valP,"",6,4) + " - " +
    // 	njm::toString(valM,"",6,4) + " -> " +
    // 	njm::toString(mu,"",6,4) + " , " + njm::toString(cm,"",6,4) +
    // 	" || " + njm::toString(par,", ","") << "\r" << std::flush;


    ++iter;
      
    mu = tp.A/std::pow(tp.B + iter,tp.ell);
    cm = tp.C/std::pow(iter,tp.t);

    
    if(mu < tp.muMin){
      converged = 1;
    }
    
  }

  agent.tp.putPar(par); // assign optimized par to the agent
}



template <class S, class A, class M, class MP>
void M1SpOptim<S,A,M,MP>
::tune(const System<M,MP,M,MP> & system,
       A agent){

  std::cout << "thread "
	    << omp_get_thread_num()
	    << " is tuning!!!!!!!" << std::endl;
  System<M,MP,M,MP> s(system.sD_r,system.tD_r,system.fD,system.dD_r,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);
  s.modelEst.fitType = MLE;
  s.fD.finalT = s.sD.time + 2*s.fD.period;

  M1SpOptim<System<M,MP,M,MP>,A,M,MP> o;
  o.tp.tune = 0;
  
  TuneRunner<System<M,MP,M,MP>,A,
	     M1SpOptim<System<M,MP,M,MP>,A,M,MP> > r;
  
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
    o.tp.A=abVals.at(i).first;
    o.tp.B=abVals.at(i).second;

    val = r.run(s,agent,o,50,s.fD.finalT);
    if(val < minVal){
      bestA = o.tp.A;
      bestB = o.tp.B;

      minVal = val;
    }
  }

  tp.A = bestA;
  tp.B = bestB;

}
