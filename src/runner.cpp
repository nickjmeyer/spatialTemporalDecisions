#include "runner.hpp"
#include "anchorMan.hpp"



template <class S, class A>
double
TrainRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  double value=0;
  int r,t;
  for(r=0; r<numReps; r++){
    system.model.assignRand(system.paramGen_r,system.paramEst_r);
    system.reset();
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      system.updateStatus();
      
      system.nextPoint();
    }
    value += system.value();
  }
  return value/((double)numReps);
}
						




template class PlainRunner<System<GravityModel,GravityParam,
				  GravityModel,GravityParam>,
			   NoTrt<GravityModel,GravityParam> >;
template class PlainRunner<System<GravityModel,GravityParam,
				  GravityModel,GravityParam>,
			   ProximalAgent<GravityModel,GravityParam> >;
template class PlainRunner<System<GravityModel,GravityParam,
				  GravityModel,GravityParam>,
			   RankAgent<ToyFeatures2<GravityModel,GravityParam>,
				     GravityModel,GravityParam> >;

template class PlainRunner<System<RangeModel,RangeParam,
				  RangeModel,RangeParam>,
			   NoTrt<RangeModel,RangeParam> >;
template class PlainRunner<System<RangeModel,RangeParam,
				  RangeModel,RangeParam>,
			   ProximalAgent<RangeModel,RangeParam> >;
template class PlainRunner<System<RangeModel,RangeParam,
				  RangeModel,RangeParam>,
			   MyopicAgent<RangeModel,RangeParam> >;
template class PlainRunner<System<RangeModel,RangeParam,
				  RangeModel,RangeParam>,
			   RankAgent<ToyFeatures2<RangeModel,RangeParam>,
				     RangeModel,RangeParam> >;

template class PlainRunner<System<GravityModel,GravityParam,
				  RangeModel,RangeParam>,
			   RankAgent<ToyFeatures2<RangeModel,RangeParam>,
				     RangeModel,RangeParam> >;

template class PlainRunner<System<CaveModel,CaveParam,
				  CaveModel,CaveParam>,
			   NoTrt<CaveModel,CaveParam> >;
template class PlainRunner<System<CaveModel,CaveParam,
				  CaveModel,CaveParam>,
			   ProximalAgent<CaveModel,CaveParam> >;

template class PlainRunner<System<CaveModel,CaveParam,
				  CaveModel,CaveParam>,
			   RankAgent<ToyFeatures2<CaveModel,CaveParam>,
				     CaveModel,CaveParam> >;




template <class S, class A>
double
PlainRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  double value=0;
  int r,t;
  for(r=0; r<numReps; r++){
    if(system.modelGen.fitType == MCMC){
      system.modelGen.mcmc.samples.setRand();
      system.paramGen_r.putPar(system.modelGen.mcmc.samples.getPar());
    }
    
    system.reset();
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      system.updateStatus();
      
      system.nextPoint();
    }
    value += system.value();
  }

  return value/((double)numReps);
}



template <class S, class A>
std::pair<double,double>
PlainRunner<S,A>
::runEx(S system,
	A agent,
	const int numReps, const int numPoints){
  double value=0,valueSq=0;
  int r,t;
  for(r=0; r<numReps; r++){
    system.reset();
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      system.updateStatus();
      
      system.nextPoint();
    }
    value += system.value();
    valueSq += system.value()*system.value();
  }
  value/=(double)numReps;
  valueSq/=(double)numReps;

  double var = (numReps/(numReps-1))*(valueSq - value*value);
  return std::pair<double,double>(value,var);
}


template class VanillaRunner<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     NoTrt<GravityModel,GravityParam> >;
template class VanillaRunner<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     ProximalAgent<GravityModel,GravityParam> >;
template class VanillaRunner<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     MyopicAgent<GravityModel,GravityParam> >;



template class VanillaRunner<System<GravityModel,GravityParam,
				    RangeModel,RangeParam>,
			     NoTrt<RangeModel,RangeParam> >;
template class VanillaRunner<System<GravityModel,GravityParam,
				    RangeModel,RangeParam>,
			     ProximalAgent<RangeModel,RangeParam> >;

template class VanillaRunner<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     RankAgent<ToyFeatures2<GravityModel,
						    GravityParam>,
				       GravityModel,GravityParam> >;



template<class S, class A>
double
VanillaRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  resetRandomSeed();
  
  double value=0;
  int r,t;

  std::vector<std::vector<double> > valueAll(numReps);
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  shared(value,valueAll)					\
  firstprivate(system,agent)					\
  private(r,t)
  for(r=0; r<numReps; r++){
    system.reset();
#pragma omp critical
    {    
      valueAll.at(r).clear();
      valueAll.at(r).push_back(system.value());
    }
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }

#pragma omp critical
    {
      value += system.value();
    }
    
    system.sD.history.push_back(system.sD.status);    
    njm::toFile(njm::toString(system.sD.history,"\n","")
		,njm::sett.datExt(agent.name+
				  "_history_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));
  }
  njm::toFile(njm::toString(valueAll,"\n",""),
	      njm::sett.datExt(agent.name+
			       "_values_",".txt"));
  return value/((double)numReps);
}




template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      GravityModel,GravityParam>,
			       NoTrt<GravityModel,GravityParam> >;
template class VanillaRunnerNS<System<GravityTimeInfModel,GravityTimeInfParam,
				      GravityTimeInfModel,GravityTimeInfParam>,
			       NoTrt<GravityTimeInfModel,GravityTimeInfParam> >;

template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      GravityModel,GravityParam>,
			       RandomAgent<GravityModel,GravityParam> >;
template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      GravityModel,GravityParam>,
			       ProximalAgent<GravityModel,GravityParam> >;
template class VanillaRunnerNS<System<GravityTimeInfModel,GravityTimeInfParam,
				      GravityTimeInfModel,GravityTimeInfParam>,
			       ProximalAgent<GravityTimeInfModel,
					     GravityTimeInfParam> >;
template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      GravityModel,GravityParam>,
			       MyopicAgent<GravityModel,GravityParam> >;



template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      RangeModel,RangeParam>,
			       NoTrt<RangeModel,RangeParam> >;
template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      RangeModel,RangeParam>,
			       ProximalAgent<RangeModel,RangeParam> >;

template class VanillaRunnerNS<System<GravityModel,GravityParam,
				      GravityModel,GravityParam>,
			       RankAgent<ToyFeatures2<GravityModel,
						      GravityParam>,
					 GravityModel,GravityParam> >;



template<class S, class A>
double
VanillaRunnerNS<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  resetRandomSeed();
  
  double value=0;
  int r,t;

#pragma omp parallel for num_threads(omp_get_max_threads())	\
  shared(value)							\
  firstprivate(system,agent)					\
  private(r,t)
  for(r=0; r<numReps; r++){
    system.reset();
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      
      system.updateStatus();
      
      system.nextPoint();

    }

#pragma omp critical
    {
      value += system.value();
    }

  }
  return value/((double)numReps);
}






template class FitOnlyRunner<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     MyopicAgent<GravityModel,GravityParam> >;
template class FitOnlyRunner<System<CaveModel,CaveParam,
				    CaveModel,CaveParam>,
			     MyopicAgent<CaveModel,CaveParam> >;
template class FitOnlyRunner<System<GravityModel,GravityParam,
				    GravityModel,GravityParam>,
			     RankAgent<ToyFeatures2<GravityModel,
						    GravityParam>,
				       GravityModel,GravityParam> >;


template class FitOnlyRunner<System<GravityModel,GravityParam,
				    RangeModel,RangeParam>,
			     MyopicAgent<RangeModel,RangeParam> >;
template class FitOnlyRunner<System<GravityModel,GravityParam,
				    CaveModel,CaveParam>,
			     MyopicAgent<CaveModel,CaveParam> >;


template <class S, class A>
double
FitOnlyRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  resetRandomSeed();
  
  double value=0;
  int r,t;

  std::vector<std::vector<double> > valueAll(numReps);
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  shared(value,valueAll)					\
  firstprivate(system,agent)					\
  private(r,t)
  for(r=0; r<numReps; r++){
    system.reset();
#pragma omp critical
    {    
      valueAll.at(r).clear();
      valueAll.at(r).push_back(system.value());
    }
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart &&
	 (((t-system.fD.trtStart) % system.fD.period) ==0))
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    system.paramEst);
      
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }

#pragma omp critical
    {
      value += system.value();
    }

    system.sD.history.push_back(system.sD.status);
    njm::toFile(njm::toString(system.sD.history,"\n","")
		,njm::sett.datExt(agent.name+
				  "_history_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));
  }
  njm::toFile(njm::toString(valueAll,"\n",""),
	      njm::sett.datExt(agent.name+
			       "_values_",".txt"));
  return value/((double)numReps);
}


template class
OptimRunner<System<GravityModel,GravityParam,
		   GravityModel,GravityParam>,
	    RankAgent<ToyFeatures2<GravityModel,GravityParam>,
		      GravityModel,GravityParam>,
	    M1SpOptim<System<GravityModel,GravityParam,
			     GravityModel,GravityParam>,
		      RankAgent<ToyFeatures2<GravityModel,
					     GravityParam>,
				GravityModel,GravityParam>,
		      GravityModel,GravityParam> >;


// Range model (misspecified)



template <class S, class A, class Optim>
double
OptimRunner<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints){

  resetRandomSeed();
  
  int tick,tickR,tock,tockR,done=0;
  tick = std::time(NULL);
  double hours;
  
  double value=0;
  int r,t;
  std::vector<std::vector<double> > valueAll(numReps);
  std::vector<std::vector<double> > weights;

  std::vector<int> times(numReps);
  
  // int threads = (omp_get_max_threads() < 16 ? 1 : omp_get_max_threads());
  int threads = omp_get_max_threads();
  
#pragma omp parallel for num_threads(threads)	\
  shared(value,valueAll,tock,tick)		\
  firstprivate(system,agent,optim,weights)	\
  private(r,t,tockR,tickR)
  for(r=0; r<numReps; r++){
    // record time for each replication
    tickR=std::time(NULL);

    system.reset();
    agent.reset();
    optim.reset();
    
#pragma omp critical
    {    
      valueAll.at(r).clear();
      valueAll.at(r).push_back(system.value());
    }
    weights.clear();

    // begin rep r
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart){
	if(t==system.fD.trtStart)
	  system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			      system.paramEst);
	else
	  system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			      system.paramEst,system.paramEst);

	
	optim.optim(system,agent);
	
	weights.push_back(agent.tp.getPar());
	
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      }
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }
    // end rep r

#pragma omp critical
    {
      value += system.value();
    }

    // store time for iteration
    tockR = std::time(NULL);
#pragma omp critical
    {
      times.at(r) = tockR - tickR;
    }

    // write history to file
    system.sD.history.push_back(system.sD.status);
    njm::toFile(njm::toString(system.sD.history,"\n","")
		,njm::sett.datExt(agent.name+"_"+optim.name+
				  "_history_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));
    // write weights
    njm::toFile(njm::toString(weights,"\n","")
		,njm::sett.datExt(agent.name+"_"+optim.name+
				  "_weights_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));


#pragma omp critical
    {
      done++;
      tock = std::time(NULL);
      hours = ((double)(tock-tick))/((double)done);
      hours /= 3600.0;

      njm::toFile("Completed " + njm::toString(done,"",6,0) +
		  " out of " + njm::toString(numReps,"",6,0) +
		  " in " + njm::toString(hours,"",8,4) + " hours" +
		  " with value " + njm::toString(value/((double)done),"",6,4) +
		  "\n",
		  njm::sett.datExt(agent.name+"_"+optim.name+"_status_",
				   ".txt"));
    }
    
  }

  njm::toFile(njm::toString(valueAll,"\n",""),
	      njm::sett.datExt(agent.name+"_"+optim.name+
			       "_values_",".txt"));

  njm::toFile(njm::toString(times,"\n",""),
	      njm::sett.datExt(agent.name+"_"+optim.name+
			       "_times_",".txt"));


  return value/((double)numReps);
}



template class
OptimRunnerNS<System<GravityModel,GravityParam,
		     GravityModel,GravityParam>,
	      RankAgent<ToyFeatures2<GravityModel,GravityParam>,
			GravityModel,GravityParam>,
	      M1SpOptim<System<GravityModel,GravityParam,
			       GravityModel,GravityParam>,
			RankAgent<ToyFeatures2<GravityModel,
					       GravityParam>,
				  GravityModel,GravityParam>,
			GravityModel,GravityParam> >;

template class
OptimRunnerNS<System<RangeModel,RangeParam,
		     RangeModel,RangeParam>,
	      RankAgent<ToyFeatures2<RangeModel,RangeParam>,
			RangeModel,RangeParam>,
	      M1SpOptim<System<RangeModel,RangeParam,
			       RangeModel,RangeParam>,
			RankAgent<ToyFeatures2<RangeModel,
					       RangeParam>,
				  RangeModel,RangeParam>,
			RangeModel,RangeParam> >;

template class
OptimRunnerNS<System<CaveModel,CaveParam,
		     CaveModel,CaveParam>,
	      RankAgent<ToyFeatures2<CaveModel,CaveParam>,
			CaveModel,CaveParam>,
	      M1SpOptim<System<CaveModel,CaveParam,
			       CaveModel,CaveParam>,
			RankAgent<ToyFeatures2<CaveModel,
					       CaveParam>,
				  CaveModel,CaveParam>,
			CaveModel,CaveParam> >;



template <class S, class A, class Optim>
double
OptimRunnerNS<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints){

  double value=0;
  int r,t;
  for(r=0; r<numReps; r++){
    system.reset();
    agent.tp.weights.ones();
    
    if(system.modelGen.fitType == MCMC){
      system.modelGen.mcmc.samples.setRand();
      system.paramGen.putPar(system.modelGen.mcmc.samples.getPar());
    }

    // begin rep r
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart &&
	 (((t-system.fD.trtStart) % system.fD.period) == 0)){
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    system.paramEst);
	optim.optim(system,agent);
      }
      
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      
      system.updateStatus();
      
      system.nextPoint();

    }
    // end rep r

    value += system.value();
  }

  return value/((double)numReps);
}







template class
TuneRunner<System<GravityModel,GravityParam,
		  GravityModel,GravityParam>,
	   RankAgent<ToyFeatures2<GravityModel,GravityParam>,
		     GravityModel,GravityParam>,
	   M1SpOptim<System<GravityModel,GravityParam,
			    GravityModel,GravityParam>,
		     RankAgent<ToyFeatures2<GravityModel,
					    GravityParam>,
			       GravityModel,GravityParam>,
		     GravityModel,GravityParam> >;

template class
TuneRunner<System<RangeModel,RangeParam,
		  RangeModel,RangeParam>,
	   RankAgent<ToyFeatures2<RangeModel,RangeParam>,
		     RangeModel,RangeParam>,
	   M1SpOptim<System<RangeModel,RangeParam,
			    RangeModel,RangeParam>,
		     RankAgent<ToyFeatures2<RangeModel,
					    RangeParam>,
			       RangeModel,RangeParam>,
		     RangeModel,RangeParam> >;

template class
TuneRunner<System<CaveModel,CaveParam,
		  CaveModel,CaveParam>,
	   RankAgent<ToyFeatures2<CaveModel,CaveParam>,
		     CaveModel,CaveParam>,
	   M1SpOptim<System<CaveModel,CaveParam,
			    CaveModel,CaveParam>,
		     RankAgent<ToyFeatures2<CaveModel,
					    CaveParam>,
			       CaveModel,CaveParam>,
		     CaveModel,CaveParam> >;




template <class S, class A, class Optim>
double
TuneRunner<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints){

  double value=0;
  int r,t;
  for(r=0; r<numReps; r++){
    if(system.modelGen.fitType == MCMC){
      system.modelGen.mcmc.samples.setRand();
      system.paramGen_r.putPar(system.modelGen.mcmc.samples.getPar());
    }
    system.reset();
    agent.tp.weights.ones();


    // begin rep r
    for(t=system.sD.time; t<numPoints; t++){
      
      if(t>=system.fD.trtStart){
	if(t==system.fD.trtStart)
	  system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			      system.paramEst);
	else
	  system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			      system.paramEst,system.paramEst);
	  
	optim.optim(system,agent);

	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      }
      
      system.updateStatus();
      
      system.nextPoint();

    }
    // end rep r

    value += system.value();
  }

  return value/((double)numReps);
}






template <class S, class A, class Optim>
double
TestRunner<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints){
  double value=0;
  int r,t;
  std::vector<std::vector<double> > valueAll(numReps);
  std::vector<std::vector<double> > weights;
  int threads = (omp_get_max_threads() < 16 ? 1 : omp_get_max_threads());
#pragma omp parallel for num_threads(threads)	\
  shared(value,valueAll)			\
  firstprivate(system,agent,optim,weights)	\
  private(r,t)
  for(r=0; r<numReps; r++){
    system.reset();
#pragma omp critical
    {    
      valueAll.at(r).clear();
      valueAll.at(r).push_back(system.value());
    }
    weights.clear();
    for(t=system.sD.time; t<numPoints; t++){
      if(t==system.fD.trtStart){
	optim.optim(system,agent);
	weights.push_back(agent.tp.getPar());	
      }
      
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }

#pragma omp critical
    {
      value += system.value();
    }
    
    system.sD.history.push_back(system.sD.status);
    njm::toFile(njm::toString(system.sD.history,"\n","")
		,njm::sett.datExt(agent.name+"_"+optim.name+
				  "_history_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));
    njm::toFile(njm::toString(weights,"\n","")
		,njm::sett.datExt(agent.name+"_"+optim.name+
				  "_weights_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));
  }
  njm::toFile(njm::toString(valueAll,"\n",""),
	      njm::sett.datExt(agent.name+"_"+optim.name+
			       "_values_",".txt"));
  return value/((double)numReps);
}





template <class S, class A>
double
TimerRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  resetRandomSeed();

  std::chrono::milliseconds fitTime,trtTime,simTime,diff;
  fitTime=std::chrono::milliseconds::zero();
  trtTime=std::chrono::milliseconds::zero();
  simTime=std::chrono::milliseconds::zero();

  int i,pointsToTime=numPoints - system.fD.trtStart;
  std::vector<std::chrono::milliseconds> trtVec(pointsToTime),
    simVec(pointsToTime);
  for(i=0; i<pointsToTime; i++){
    trtVec.at(i) = std::chrono::milliseconds::zero();
    simVec.at(i) = std::chrono::milliseconds::zero();
  }
  
  std::chrono::time_point< std::chrono::high_resolution_clock > tick,tock;
  int r,t;

  for(r=0; r<numReps; r++){
    system.reset();

    for(t=system.sD.time; t<numPoints; t++){
      // fit
      tick=std::chrono::high_resolution_clock::now();
      if(t>=system.fD.trtStart &&
	 (((t-system.fD.trtStart) % system.fD.period) ==0))
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    system.paramEst);
      tock=std::chrono::high_resolution_clock::now();
      
      diff=std::chrono::milliseconds::zero();
      diff += std::chrono::duration_cast< std::chrono::milliseconds >
	(tock.time_since_epoch());
      diff -= std::chrono::duration_cast< std::chrono::milliseconds >
	(tick.time_since_epoch());
      fitTime += diff;
	  

      // trt
      tick=std::chrono::high_resolution_clock::now();
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      tock=std::chrono::high_resolution_clock::now();

      diff=std::chrono::milliseconds::zero();
      diff += std::chrono::duration_cast< std::chrono::milliseconds >
	(tock.time_since_epoch());
      diff -= std::chrono::duration_cast< std::chrono::milliseconds >
	(tick.time_since_epoch());
      trtTime += diff;

      for(i=t-system.fD.trtStart; i>=0 && i<pointsToTime; i++)
	trtVec.at(i) += diff;

      // sim
      tick=std::chrono::high_resolution_clock::now();
      system.updateStatus();
      
      system.nextPoint();
      tock=std::chrono::high_resolution_clock::now();

      diff=std::chrono::milliseconds::zero();
      diff += std::chrono::duration_cast< std::chrono::milliseconds >
	(tock.time_since_epoch());
      diff -= std::chrono::duration_cast< std::chrono::milliseconds >
	(tick.time_since_epoch());
      simTime += diff;

      for(i=t-system.fD.trtStart; i>=0 && i<pointsToTime; i++)
	simVec.at(i) += diff;
    }

  }

  double fitTotal,trtTotal,simTotal;
  fitTotal = (1.6666667e-5 * ((double)fitTime.count()))/((double)numReps);
  trtTotal = (1.6666667e-5 * ((double)trtTime.count()))/((double)numReps);
  simTotal = (1.6666667e-5 * ((double)simTime.count()))/((double)numReps);

  double total = fitTotal+trtTotal+simTotal;
  double hours = 1.6666667e-2 * ((double)total);
  
  std::string times=njm::toString("fit",": ",5,0) +
    njm::toString(fitTotal,"minutes\n",12,6) +
    njm::toString("trt",": ",5,0) +
    njm::toString(trtTotal,"minutes\n",12,6) +
    njm::toString("sim",": ",5,0) +
    njm::toString(simTotal,"minutes\n",12,6) +
    "\n" +
    njm::toString("tot",": ",5,0) +
    njm::toString(hours,"hours",12,6);

  times += "\n\n";
  times += njm::toString("t","",5,0);
  times += njm::toString("trt","",12,0);
  times += njm::toString("sim","",12,0);
  times += "\n";
  
  for(i=0; i<pointsToTime; i++){
    times += njm::toString(i + system.fD.trtStart,"",5,0);
    
    total = (1.6666667e-5 * ((double)trtVec.at(i).count()))/((double)numReps);
    times += njm::toString(total,"",12,6);

    total = (1.6666667e-5 * ((double)simVec.at(i).count()))/((double)numReps);
    times += njm::toString(total,"",12,6);

    times += njm::toString("minutes","",12,0);
    
    times += "\n";
  }

  njm::message(times);

  return (hours);
}
