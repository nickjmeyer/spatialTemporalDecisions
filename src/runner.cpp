#include "runner.hpp"
#include "anchorMan.hpp"



template <class S, class A>
RunStats
TrainRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  // double value=0;
  int r,t;
  RunStats rs;
  for(r=0; r<numReps; r++){
    system.model.assignRand(system.paramGen_r,system.paramEst_r);
    system.revert();
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      system.updateStatus();
      
      system.nextPoint();
    }
    rs(system.value());
    // value += system.value();
  }
  return rs;
}
						





template class PlainRunner<System<ModelTimeExpCaves,
				  ModelTimeExpCaves>,
			   RankAgent<ToyFeatures4<ModelTimeExpCaves>,
				     ModelTimeExpCaves> >;



template class PlainRunner<System<ModelTimeExpCaves,
				  ModelTimeExpCaves>,
			   RankAgent<ToyFeatures5<ModelTimeExpCaves>,
				     ModelTimeExpCaves> >;


template class PlainRunner<System<ModelRadius,
				  ModelRadius>,
			   RankAgent<ToyFeatures4<ModelRadius>,
				     ModelRadius> >;



template class PlainRunner<System<ModelRadius,
				  ModelRadius>,
			   RankAgent<ToyFeatures5<ModelRadius>,
				     ModelRadius> >;


template class PlainRunner<System<ModelDist,
				  ModelDist>,
			   RankAgent<ToyFeatures4<ModelDist>,
				     ModelDist> >;

template class PlainRunner<System<ModelDist,
				  ModelDist>,
			   RankAgent<ToyFeatures5<ModelDist>,
				     ModelDist> >;


template class PlainRunner<System<ModelCovar,
				  ModelCovar>,
			   RankAgent<ToyFeatures5<ModelCovar>,
				     ModelCovar> >;


template class PlainRunner<System<ModelDistKern,
				  ModelDistKern>,
			   RankAgent<ToyFeatures4<ModelDistKern>,
				     ModelDistKern> >;



template class PlainRunner<System<ModelDistKern,
				  ModelDistKern>,
			   RankAgent<ToyFeatures5<ModelDistKern>,
				     ModelDistKern> >;


template class PlainRunner<System<ModelGravity,
				  ModelGravity>,
			   NoTrt<ModelGravity> >;


template class PlainRunner<System<ModelTimeExpCaves,
				  ModelTimeExpCaves>,
			   RankAgent<WnsFeatures1<ModelTimeExpCaves>,
				     ModelTimeExpCaves> >;


template class PlainRunner<System<ModelGravity,
				  ModelGravity>,
			   RankAgent<WnsFeatures1<ModelGravity>,
				     ModelGravity> >;


template class PlainRunner<System<ModelRadius,
				  ModelRadius>,
			   RankAgent<WnsFeatures1<ModelRadius>,
				     ModelRadius> >;


template class PlainRunner<System<ModelDist,
				  ModelDist>,
			   RankAgent<WnsFeatures1<ModelDist>,
				     ModelDist> >;


template class PlainRunner<System<ModelDistKern,
				  ModelDistKern>,
			   RankAgent<WnsFeatures1<ModelDistKern>,
				     ModelDistKern> >;





template <class S, class A>
RunStats
PlainRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints){
  RunStats rs;
  // double value=0;
  int r,t;
  for(r=0; r<numReps; r++){
    if(system.modelGen_r.sample()){
      std::vector<double> newPar = system.modelGen_r.getPar();
      system.modelEst_r.putPar(newPar.begin());
      system.modelGen_r.setFill(system.sD,system.tD,system.fD,system.dD);
      system.modelEst_r.setFill(system.sD,system.tD,system.fD,system.dD);
    }
    
    system.revert();

    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);

      system.updateStatus();

      system.nextPoint();
    }
    
    rs(system.value());
    // value += system.value();
  }

  return rs;
}



template class VanillaRunner<System<ModelTimeExpCaves,
				    ModelTimeExpCaves>,
			     NoTrt<ModelTimeExpCaves> >;

template class VanillaRunner<System<ModelTimeExpCaves,
				    ModelTimeExpCaves>,
			     ProximalAgent<ModelTimeExpCaves> >;

template class VanillaRunner<System<ModelGravity,
				    ModelGravity>,
			     NoTrt<ModelGravity> >;

template class VanillaRunner<System<ModelGravity,
				    ModelGravity>,
			     ProximalAgent<ModelGravity> >;

template class VanillaRunner<System<ModelTimeExpCaves,
				    ModelRadius>,
			     NoTrt<ModelRadius> >;

template class VanillaRunner<System<ModelTimeExpCaves,
				    ModelRadius>,
			     ProximalAgent<ModelRadius> >;



template<class S, class A>
RunStats
VanillaRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints,
      const Starts & starts){
  // double value=0;
  int r,t;

  RunStats rs;
  std::vector<std::vector<double> > valueAll(numReps);
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  shared(valueAll,starts,rs)					\
  firstprivate(system,agent)					\
  private(r,t)
  for(r=0; r<numReps; r++){
    njm::resetSeed(r);
    system.reset(starts[r]);
    
#pragma omp critical
    {    
      valueAll.at(r).clear();
      valueAll.at(r).push_back(system.value());
    }
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }

#pragma omp critical
    {
      rs(system.value());
      // value += system.value();
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
  return rs;
}



template class VanillaRunnerNS<System<ModelTimeExpCaves,
				      ModelTimeExpCaves>,
			       NoTrt<ModelTimeExpCaves> >;

template class VanillaRunnerNS<System<ModelTimeExpCaves,
				      ModelTimeExpCaves>,
			       ProximalAgent<ModelTimeExpCaves> >;

template class VanillaRunnerNS<System<ModelTimeExpCaves,
				      ModelTimeExpCaves>,
			       MyopicAgent<ModelTimeExpCaves> >;

template class VanillaRunnerNS<System<ModelGravity,
				      ModelGravity>,
			       NoTrt<ModelGravity> >;

template class VanillaRunnerNS<System<ModelGravity,
				      ModelGravity>,
			       ProximalAgent<ModelGravity> >;

template class VanillaRunnerNS<System<ModelGravity,
				      ModelGravity>,
			       MyopicAgent<ModelGravity> >;

template class VanillaRunnerNS<System<ModelGravity,
				      ModelGravity>,
			       RankAgent<ToyFeatures4<ModelGravity>,
					 ModelGravity> >;

template class VanillaRunnerNS<System<ModelTimeExpCaves,
				      ModelTimeExpCaves>,
			       RankAgent<ToyFeatures4<ModelTimeExpCaves>,
					 ModelTimeExpCaves> >;

template class VanillaRunnerNS<System<ModelGravity,
				      ModelGravity>,
			       RankAgent<ToyFeatures5<ModelGravity>,
					 ModelGravity> >;

template class VanillaRunnerNS<System<ModelTimeExpCaves,
				      ModelTimeExpCaves>,
			       RankAgent<ToyFeatures5<ModelTimeExpCaves>,
					 ModelTimeExpCaves> >;

template class VanillaRunnerNS<System<ModelGravity,
				      ModelGravity>,
			       RankAgent<WnsFeatures1<ModelGravity>,
					 ModelGravity> >;

template class VanillaRunnerNS<System<ModelTimeExpCaves,
				      ModelTimeExpCaves>,
			       RankAgent<WnsFeatures1<ModelTimeExpCaves>,
					 ModelTimeExpCaves> >;





template<class S, class A>
RunStats
VanillaRunnerNS<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints,
      const Starts & starts){
  // double value=0;
  int r,t;

  RunStats rs;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  shared(starts,rs)						\
  firstprivate(system,agent)					\
  private(r,t)
  for(r=0; r<numReps; r++){
    njm::resetSeed(r);
    system.reset(starts[r]);
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);
      
      system.updateStatus();
      
      system.nextPoint();

    }

#pragma omp critical
    {
      rs(system.value());
      // value += system.value();
    }

  }
  return rs;
}






template class FitOnlyRunner<System<ModelTimeExpCaves,
				    ModelTimeExpCaves>,
			     MyopicAgent<ModelTimeExpCaves> >;

template class FitOnlyRunner<System<ModelTimeExpCaves,
				    ModelRadius>,
			     MyopicAgent<ModelRadius> >;

template class FitOnlyRunner<System<ModelTimeExpCaves,
				    ModelDist>,
			     MyopicAgent<ModelDist> >;

template class FitOnlyRunner<System<ModelTimeExpCaves,
				    ModelCovar>,
			     MyopicAgent<ModelCovar> >;

template class FitOnlyRunner<System<ModelTimeExpCaves,
				    ModelDistKern>,
			     MyopicAgent<ModelDistKern> >;



template <class S, class A>
RunStats
FitOnlyRunner<S,A>
::run(S system,
      A agent,
      const int numReps, const int numPoints,
      const Starts & starts){
  // double value=0;
  int r,t;

  RunStats rs;
  std::vector<std::vector<double> > valueAll(numReps);
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  shared(valueAll,starts,rs)					\
  firstprivate(system,agent)					\
  private(r,t)
  for(r=0; r<numReps; r++){
    njm::resetSeed(r);
    system.reset(starts[r]);
    
#pragma omp critical
    {    
      valueAll.at(r).clear();
      valueAll.at(r).push_back(system.value());
    }
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart){
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    t > system.fD.trtStart);

	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);
      }
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }

#pragma omp critical
    {
      rs(system.value());
      // value += system.value();
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
  return rs;
}


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelTimeExpCaves>,
	    RankAgent<ToyFeatures4<ModelTimeExpCaves>,
		      ModelTimeExpCaves>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelTimeExpCaves>,
		      RankAgent<ToyFeatures4<ModelTimeExpCaves>,
				ModelTimeExpCaves>,
		      ModelTimeExpCaves> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelRadius>,
	    RankAgent<ToyFeatures4<ModelRadius>,
		      ModelRadius>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelRadius>,
		      RankAgent<ToyFeatures4<ModelRadius>,
				ModelRadius>,
		      ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDist>,
	    RankAgent<ToyFeatures4<ModelDist>,
		      ModelDist>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelDist>,
		      RankAgent<ToyFeatures4<ModelDist>,
				ModelDist>,
		      ModelDist> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDistKern>,
	    RankAgent<ToyFeatures4<ModelDistKern>,
		      ModelDistKern>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelDistKern>,
		      RankAgent<ToyFeatures4<ModelDistKern>,
				ModelDistKern>,
		      ModelDistKern> >;




template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelTimeExpCaves>,
	    RankAgent<ToyFeatures5<ModelTimeExpCaves>,
		      ModelTimeExpCaves>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelTimeExpCaves>,
		      RankAgent<ToyFeatures5<ModelTimeExpCaves>,
				ModelTimeExpCaves>,
		      ModelTimeExpCaves> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelRadius>,
	    RankAgent<ToyFeatures5<ModelRadius>,
		      ModelRadius>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelRadius>,
		      RankAgent<ToyFeatures5<ModelRadius>,
				ModelRadius>,
		      ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDist>,
	    RankAgent<ToyFeatures5<ModelDist>,
		      ModelDist>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelDist>,
		      RankAgent<ToyFeatures5<ModelDist>,
				ModelDist>,
		      ModelDist> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelCovar>,
	    RankAgent<ToyFeatures5<ModelCovar>,
		      ModelCovar>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelCovar>,
		      RankAgent<ToyFeatures5<ModelCovar>,
				ModelCovar>,
		      ModelCovar> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDistKern>,
	    RankAgent<ToyFeatures5<ModelDistKern>,
		      ModelDistKern>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelDistKern>,
		      RankAgent<ToyFeatures5<ModelDistKern>,
				ModelDistKern>,
		      ModelDistKern> >;






template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelTimeExpCaves>,
	    RankAgent<WnsFeatures1<ModelTimeExpCaves>,
		      ModelTimeExpCaves>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelTimeExpCaves>,
		      RankAgent<WnsFeatures1<ModelTimeExpCaves>,
				ModelTimeExpCaves>,
		      ModelTimeExpCaves> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelRadius>,
	    RankAgent<WnsFeatures1<ModelRadius>,
		      ModelRadius>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelRadius>,
		      RankAgent<WnsFeatures1<ModelRadius>,
				ModelRadius>,
		      ModelRadius> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDist>,
	    RankAgent<WnsFeatures1<ModelDist>,
		      ModelDist>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelDist>,
		      RankAgent<WnsFeatures1<ModelDist>,
				ModelDist>,
		      ModelDist> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDistKern>,
	    RankAgent<WnsFeatures1<ModelDistKern>,
		      ModelDistKern>,
	    M1SpOptim<System<ModelTimeExpCaves,
			     ModelDistKern>,
		      RankAgent<WnsFeatures1<ModelDistKern>,
				ModelDistKern>,
		      ModelDistKern> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelTimeExpCaves>,
	    OsspAgent<ModelTimeExpCaves>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelTimeExpCaves>,
			OsspAgent<ModelTimeExpCaves>,
			ToyFeatures4<ModelTimeExpCaves>,
			ModelTimeExpCaves> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelRadius>,
	    OsspAgent<ModelRadius>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelRadius>,
			OsspAgent<ModelRadius>,
			ToyFeatures4<ModelRadius>,
			ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDist>,
	    OsspAgent<ModelDist>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelDist>,
			OsspAgent<ModelDist>,
			ToyFeatures4<ModelDist>,
			ModelDist> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelTimeExpCaves>,
	    OsspAgent<ModelTimeExpCaves>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelTimeExpCaves>,
			OsspAgent<ModelTimeExpCaves>,
			ToyFeatures5<ModelTimeExpCaves>,
			ModelTimeExpCaves> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelRadius>,
	    OsspAgent<ModelRadius>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelRadius>,
			OsspAgent<ModelRadius>,
			ToyFeatures5<ModelRadius>,
			ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDist>,
	    OsspAgent<ModelDist>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelDist>,
			OsspAgent<ModelDist>,
			ToyFeatures5<ModelDist>,
			ModelDist> >;


template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelCovar>,
	    OsspAgent<ModelCovar>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelCovar>,
			OsspAgent<ModelCovar>,
			ToyFeatures5<ModelCovar>,
			ModelCovar> >;



template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelTimeExpCaves>,
	    OsspAgent<ModelTimeExpCaves>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelTimeExpCaves>,
			OsspAgent<ModelTimeExpCaves>,
			WnsFeatures1<ModelTimeExpCaves>,
			ModelTimeExpCaves> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelRadius>,
	    OsspAgent<ModelRadius>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelRadius>,
			OsspAgent<ModelRadius>,
			WnsFeatures1<ModelRadius>,
			ModelRadius> >;

template class
OptimRunner<System<ModelTimeExpCaves,
		   ModelDist>,
	    OsspAgent<ModelDist>,
	    M1OsspOptim<System<ModelTimeExpCaves,
			       ModelDist>,
			OsspAgent<ModelDist>,
			WnsFeatures1<ModelDist>,
			ModelDist> >;



template <class S, class A, class Optim>
RunStats
OptimRunner<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints,
      const Starts & starts){
  int tick,tickR,tock,tockR,done=0;
  tick = std::time(NULL);
  double hours;

  RunStats rs;
  // double value=0;
  int r,t;
  std::vector<std::vector<double> > valueAll(numReps);
  std::vector<std::vector<double> > weights;

  std::vector<int> times(numReps);
  
  // int threads = (omp_get_max_threads() < 16 ? 1 : omp_get_max_threads());
  int threads = omp_get_max_threads();
  
#pragma omp parallel for num_threads(threads)	\
  shared(valueAll,tock,tick,starts,rs)		\
  firstprivate(system,agent,optim,weights)	\
  private(r,t,tockR,tickR)
  for(r=0; r<numReps; r++){
    njm::resetSeed(r);
    // record time for each replication
    tickR=std::time(NULL);

    system.reset(starts[r]);
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
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    t > system.fD.trtStart);
	
	optim.optim(system,agent);
	
	weights.push_back(agent.tp.getPar());
	
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);
      }
      
      system.updateStatus();
      
      system.nextPoint();

#pragma omp critical
      {
	valueAll.at(r).push_back(system.value());
      }
    }
    // end rep r...time to write results to disk

#pragma omp critical
    {
      rs(system.value());
      // value += system.value();
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

    // write optim parameters to file
    njm::toFile(njm::toString(optim.tp.getPar()," ","\n")
		,njm::sett.datExt(agent.name+"_"+optim.name+
				  "_tunePar_"+
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
		  " with value " + njm::toString(rs.smean(),"",6,4) +
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


  return rs;
}



template <class S, class A, class Optim>
RunStats
OptimRunnerNS<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints,
      const Starts & starts){

  RunStats rs;
  // double value=0;
  int r,t;
  
  int threads = omp_get_max_threads();
  
#pragma omp parallel for num_threads(threads)	\
  shared(starts,rs)				\
  firstprivate(system,agent,optim)		\
  private(r,t)
  for(r=0; r<numReps; r++){
    system.reset(starts[r]);
    agent.reset();
    optim.reset();
    
    // begin rep r
    for(t=system.sD.time; t<numPoints; t++){
      if(t>=system.fD.trtStart){
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    t > system.fD.trtStart);

	optim.optim(system,agent);
	
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);
      }
      
      system.updateStatus();
      
      system.nextPoint();

    }
    // end rep r

#pragma omp critical
    {
      rs(system.value());
      // value += system.value();
    }

  }

  return rs;
}








template class
TuneRunner<System<ModelTimeExpCaves,
		  ModelTimeExpCaves>,
	   RankAgent<ToyFeatures4<ModelTimeExpCaves>,
		     ModelTimeExpCaves>,
	   M1SpOptim<System<ModelTimeExpCaves,
			    ModelTimeExpCaves>,
		     RankAgent<ToyFeatures4<ModelTimeExpCaves>,
			       ModelTimeExpCaves>,
		     ModelTimeExpCaves> >;


template class
TuneRunner<System<ModelRadius,
		  ModelRadius>,
	   RankAgent<ToyFeatures4<ModelRadius>,
		     ModelRadius>,
	   M1SpOptim<System<ModelRadius,
			    ModelRadius>,
		     RankAgent<ToyFeatures4<ModelRadius>,
			       ModelRadius>,
		     ModelRadius> >;


template class
TuneRunner<System<ModelDist,
		  ModelDist>,
	   RankAgent<ToyFeatures4<ModelDist>,
		     ModelDist>,
	   M1SpOptim<System<ModelDist,
			    ModelDist>,
		     RankAgent<ToyFeatures4<ModelDist>,
			       ModelDist>,
		     ModelDist> >;


template class
TuneRunner<System<ModelDistKern,
		  ModelDistKern>,
	   RankAgent<ToyFeatures4<ModelDistKern>,
		     ModelDistKern>,
	   M1SpOptim<System<ModelDistKern,
			    ModelDistKern>,
		     RankAgent<ToyFeatures4<ModelDistKern>,
			       ModelDistKern>,
		     ModelDistKern> >;


template class
TuneRunner<System<ModelTimeExpCaves,
		  ModelTimeExpCaves>,
	   RankAgent<ToyFeatures5<ModelTimeExpCaves>,
		     ModelTimeExpCaves>,
	   M1SpOptim<System<ModelTimeExpCaves,
			    ModelTimeExpCaves>,
		     RankAgent<ToyFeatures5<ModelTimeExpCaves>,
			       ModelTimeExpCaves>,
		     ModelTimeExpCaves> >;


template class
TuneRunner<System<ModelRadius,
		  ModelRadius>,
	   RankAgent<ToyFeatures5<ModelRadius>,
		     ModelRadius>,
	   M1SpOptim<System<ModelRadius,
			    ModelRadius>,
		     RankAgent<ToyFeatures5<ModelRadius>,
			       ModelRadius>,
		     ModelRadius> >;


template class
TuneRunner<System<ModelDist,
		  ModelDist>,
	   RankAgent<ToyFeatures5<ModelDist>,
		     ModelDist>,
	   M1SpOptim<System<ModelDist,
			    ModelDist>,
		     RankAgent<ToyFeatures5<ModelDist>,
			       ModelDist>,
		     ModelDist> >;


template class
TuneRunner<System<ModelCovar,
		  ModelCovar>,
	   RankAgent<ToyFeatures5<ModelCovar>,
		     ModelCovar>,
	   M1SpOptim<System<ModelCovar,
			    ModelCovar>,
		     RankAgent<ToyFeatures5<ModelCovar>,
			       ModelCovar>,
		     ModelCovar> >;



template class
TuneRunner<System<ModelDistKern,
		  ModelDistKern>,
	   RankAgent<ToyFeatures5<ModelDistKern>,
		     ModelDistKern>,
	   M1SpOptim<System<ModelDistKern,
			    ModelDistKern>,
		     RankAgent<ToyFeatures5<ModelDistKern>,
			       ModelDistKern>,
		     ModelDistKern> >;



template class
TuneRunner<System<ModelTimeExpCaves,
		  ModelTimeExpCaves>,
	   RankAgent<WnsFeatures1<ModelTimeExpCaves>,
		     ModelTimeExpCaves>,
	   M1SpOptim<System<ModelTimeExpCaves,
			    ModelTimeExpCaves>,
		     RankAgent<WnsFeatures1<ModelTimeExpCaves>,
			       ModelTimeExpCaves>,
		     ModelTimeExpCaves> >;


template class
TuneRunner<System<ModelRadius,
		  ModelRadius>,
	   RankAgent<WnsFeatures1<ModelRadius>,
		     ModelRadius>,
	   M1SpOptim<System<ModelRadius,
			    ModelRadius>,
		     RankAgent<WnsFeatures1<ModelRadius>,
			       ModelRadius>,
		     ModelRadius> >;


template class
TuneRunner<System<ModelDist,
		  ModelDist>,
	   RankAgent<WnsFeatures1<ModelDist>,
		     ModelDist>,
	   M1SpOptim<System<ModelDist,
			    ModelDist>,
		     RankAgent<WnsFeatures1<ModelDist>,
			       ModelDist>,
		     ModelDist> >;




template class
TuneRunner<System<ModelDistKern,
		  ModelDistKern>,
	   RankAgent<WnsFeatures1<ModelDistKern>,
		     ModelDistKern>,
	   M1SpOptim<System<ModelDistKern,
			    ModelDistKern>,
		     RankAgent<WnsFeatures1<ModelDistKern>,
			       ModelDistKern>,
		     ModelDistKern> >;


template class
TuneRunner<System<ModelGravity,
		  ModelGravity>,
	   RankAgent<WnsFeatures1<ModelGravity>,
		     ModelGravity>,
	   M1SpOptim<System<ModelGravity,
			    ModelGravity>,
		     RankAgent<WnsFeatures1<ModelGravity>,
			       ModelGravity>,
		     ModelGravity> >;






template <class S, class A, class Optim>
RunStats
TuneRunner<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints){

  RunStats rs;
  // double value=0;
  int r,t;
  for(r=0; r<numReps; r++){
    if(system.modelGen_r.sample()){
      std::vector<double> newPar = system.modelGen_r.getPar();
      system.modelEst_r.putPar(newPar.begin());
    }

    system.revert();
    agent.reset();

    // begin rep r
    for(t=system.sD.time; t<numPoints; t++){
      
      if(t>=system.fD.trtStart){
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    t > system.fD.trtStart);
	  
	optim.optim(system,agent);

	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst);
      }
      
      system.updateStatus();
      
      system.nextPoint();

    }
    // end rep r

    rs(system.value());
    // value += system.value();
  }

  return rs;
}






template <class S, class A, class Optim>
RunStats
TestRunner<S,A,Optim>
::run(S system,
      A agent,
      Optim optim,
      const int numReps, const int numPoints){

  RunStats rs;
  // double value=0;
  int r,t;
  std::vector<std::vector<double> > valueAll(numReps);
  std::vector<std::vector<double> > weights;
  int threads = (omp_get_max_threads() < 16 ? 1 : omp_get_max_threads());
#pragma omp parallel for num_threads(threads)	\
  shared(valueAll,rs)				\
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
      rs(system.value());
      // value += system.value();
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
  return rs;
}
