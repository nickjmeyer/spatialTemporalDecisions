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
						



template class PlainRunner<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   RankAgent<ToyFeatures0<ModelTimeExpCavesGDist>,
				     ModelTimeExpCavesGDist> >;


template class PlainRunner<System<ModelTimeExpCavesEDist,
				  ModelTimeExpCavesEDist>,
			   RankAgent<ToyFeatures0<ModelTimeExpCavesEDist>,
				     ModelTimeExpCavesEDist> >;


template class PlainRunner<System<ModelGDist,
				  ModelGDist>,
			   RankAgent<ToyFeatures0<ModelGDist>,
				     ModelGDist> >;


template class PlainRunner<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   RankAgent<ToyFeatures4<ModelTimeExpCavesGDist>,
				     ModelTimeExpCavesGDist> >;



template class PlainRunner<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
				     ModelTimeExpCavesGDist> >;


template class PlainRunner<System<ModelTimeGDistTrendPow,
				  ModelTimeGDistTrendPow>,
			   RankAgent<ToyFeatures5<ModelTimeGDistTrendPow>,
				     ModelTimeGDistTrendPow> >;


template class
PlainRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    RankAgent<ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;


template class
PlainRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    RankAgent<ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;


template class
PlainRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    RankAgent<ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;


template class PlainRunner<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesEDist>,
			   RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
				     ModelTimeExpCavesEDist> >;


template class PlainRunner<System<ModelTimeExpCavesEDist,
				  ModelTimeExpCavesEDist>,
			   RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
				     ModelTimeExpCavesEDist> >;


template class PlainRunner<System<ModelRadius,
				  ModelRadius>,
			   RankAgent<ToyFeatures4<ModelRadius>,
				     ModelRadius> >;



template class PlainRunner<System<ModelRadius,
				  ModelRadius>,
			   RankAgent<ToyFeatures5<ModelRadius>,
				     ModelRadius> >;


template class PlainRunner<System<ModelGDist,
				  ModelGDist>,
			   RankAgent<ToyFeatures4<ModelGDist>,
				     ModelGDist> >;

template class PlainRunner<System<ModelGDist,
				  ModelGDist>,
			   RankAgent<ToyFeatures5<ModelGDist>,
				     ModelGDist> >;


template class PlainRunner<System<ModelCovar,
				  ModelCovar>,
			   RankAgent<ToyFeatures5<ModelCovar>,
				     ModelCovar> >;


template class PlainRunner<System<ModelGDistKern,
				  ModelGDistKern>,
			   RankAgent<ToyFeatures4<ModelGDistKern>,
				     ModelGDistKern> >;



template class PlainRunner<System<ModelGDistKern,
				  ModelGDistKern>,
			   RankAgent<ToyFeatures5<ModelGDistKern>,
				     ModelGDistKern> >;


template class PlainRunner<System<ModelGravityGDist,
				  ModelGravityGDist>,
			   NoTrt<ModelGravityGDist> >;


template class PlainRunner<System<ModelTimeExpCavesGDist,
				  ModelTimeExpCavesGDist>,
			   RankAgent<WnsFeatures1<ModelTimeExpCavesGDist>,
				     ModelTimeExpCavesGDist> >;


template class PlainRunner<System<ModelTimeExpCavesEDist,
				  ModelTimeExpCavesEDist>,
			   RankAgent<WnsFeatures1<ModelTimeExpCavesEDist>,
				     ModelTimeExpCavesEDist> >;


template class PlainRunner<System<ModelGravityGDist,
				  ModelGravityGDist>,
			   RankAgent<WnsFeatures1<ModelGravityGDist>,
				     ModelGravityGDist> >;


template class PlainRunner<System<ModelGravityEDist,
				  ModelGravityEDist>,
			   RankAgent<WnsFeatures1<ModelGravityEDist>,
				     ModelGravityEDist> >;


template class PlainRunner<System<ModelRadius,
				  ModelRadius>,
			   RankAgent<WnsFeatures1<ModelRadius>,
				     ModelRadius> >;


template class PlainRunner<System<ModelGDist,
				  ModelGDist>,
			   RankAgent<WnsFeatures1<ModelGDist>,
				     ModelGDist> >;


template class PlainRunner<System<ModelGDistKern,
				  ModelGDistKern>,
			   RankAgent<WnsFeatures1<ModelGDistKern>,
				     ModelGDistKern> >;





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



template class VanillaRunner<System<ModelTimeExpCavesGDist,
				    ModelTimeExpCavesGDist>,
			     NoTrt<ModelTimeExpCavesGDist> >;

template class VanillaRunner<System<ModelTimeGDistTrendPow,
				    ModelTimeGDistTrendPow>,
			     NoTrt<ModelTimeGDistTrendPow> >;

template class VanillaRunner<System<ModelTimeExpCavesGDistTrendPowCon,
				    ModelTimeExpCavesGDistTrendPowCon>,
			     NoTrt<ModelTimeExpCavesGDistTrendPowCon> >;

template class VanillaRunner<System<ModelTimeExpCavesGDist,
				    ModelTimeExpCavesGDist>,
			     ProximalGDistAgent<ModelTimeExpCavesGDist> >;

template class VanillaRunner<System<ModelTimeGDistTrendPow,
				    ModelTimeGDistTrendPow>,
			     ProximalGDistAgent<ModelTimeGDistTrendPow> >;

template class
VanillaRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		     ModelTimeExpCavesGDistTrendPowCon>,
	      ProximalGDistAgent<ModelTimeExpCavesGDistTrendPowCon> >;

template class
VanillaRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		     ModelTimeExpCavesEDist>,
	      ProximalEDistAgent<ModelTimeExpCavesEDist> >;

template class
VanillaRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		     ModelTimeExpCavesGDist>,
	      ProximalGDistAgent<ModelTimeExpCavesGDist> >;

template class VanillaRunner<System<ModelTimeExpCavesGDist,
				    ModelTimeExpCavesEDist>,
			     ProximalEDistAgent<ModelTimeExpCavesEDist> >;

template class VanillaRunner<System<ModelGravityGDist,
				    ModelGravityGDist>,
			     NoTrt<ModelGravityGDist> >;

template class VanillaRunner<System<ModelGravityGDist,
				    ModelGravityGDist>,
			     ProximalGDistAgent<ModelGravityGDist> >;

template class VanillaRunner<System<ModelTimeExpCavesGDist,
				    ModelRadius>,
			     NoTrt<ModelRadius> >;

template class VanillaRunner<System<ModelTimeExpCavesGDist,
				    ModelRadius>,
			     ProximalGDistAgent<ModelRadius> >;



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
		,njm::sett.datExt(agent.name +
				  "_history_"+
				  njm::toString(r,"",0,0)
				  +"_",".txt"));
  }
  njm::toFile(njm::toString(valueAll,"\n",""),
	      njm::sett.datExt(agent.name+
			       "_values_",".txt"));
  return rs;
}



template class VanillaRunnerNS<System<ModelTimeGDistTrendPow,
				      ModelTimeGDistTrendPow>,
			       NoTrt<ModelTimeGDistTrendPow> >;

template class VanillaRunnerNS<System<ModelTimeGDistTrendPow,
				      ModelTimeGDistTrendPow>,
			       ProximalGDistAgent<ModelTimeGDistTrendPow> >;

template class VanillaRunnerNS<System<ModelTimeGDistTrendPow,
				      ModelTimeGDistTrendPow>,
			       MyopicAgent<ModelTimeGDistTrendPow> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDistTrendPowCon,
				      ModelTimeExpCavesGDistTrendPowCon>,
			       NoTrt<ModelTimeExpCavesGDistTrendPowCon> >;

template class
VanillaRunnerNS<System<ModelTimeExpCavesGDistTrendPowCon,
		       ModelTimeExpCavesGDistTrendPowCon>,
		ProximalGDistAgent<ModelTimeExpCavesGDistTrendPowCon> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDistTrendPowCon,
				      ModelTimeExpCavesGDistTrendPowCon>,
			       MyopicAgent<ModelTimeExpCavesGDistTrendPowCon> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDist,
				      ModelTimeExpCavesGDist>,
			       NoTrt<ModelTimeExpCavesGDist> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDist,
				      ModelTimeExpCavesGDist>,
			       ProximalGDistAgent<ModelTimeExpCavesGDist> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDist,
				      ModelTimeExpCavesGDist>,
			       MyopicAgent<ModelTimeExpCavesGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
				      ModelGravityGDist>,
			       NoTrt<ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
				      ModelGravityGDist>,
			       ProximalGDistAgent<ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
				      ModelGravityGDist>,
			       MyopicAgent<ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
				      ModelGravityGDist>,
			       RankAgent<ToyFeatures4<ModelGravityGDist>,
					 ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelTimeGDistTrendPow,
				      ModelTimeGDistTrendPow>,
			       RankAgent<ToyFeatures4<ModelTimeGDistTrendPow>,
					 ModelTimeGDistTrendPow> >;

template class
VanillaRunnerNS<System<ModelTimeExpCavesGDistTrendPowCon,
		       ModelTimeExpCavesGDistTrendPowCon>,
		RankAgent<ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
			  ModelTimeExpCavesGDistTrendPowCon> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
				      ModelGravityGDist>,
			       RankAgent<ToyFeatures5<ModelGravityGDist>,
					 ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDist,
				      ModelTimeExpCavesGDist>,
			       RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
					 ModelTimeExpCavesGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
				      ModelGravityGDist>,
			       RankAgent<WnsFeatures1<ModelGravityGDist>,
					 ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDist,
				      ModelTimeExpCavesGDist>,
			       RankAgent<WnsFeatures1<ModelTimeExpCavesGDist>,
					 ModelTimeExpCavesGDist> >;

template class VanillaRunnerNS<System<ModelTimeExpCavesGDist,
				      ModelTimeExpCavesGDist>,
			       RankAgent<ToyFeatures7<ModelTimeExpCavesGDist>,
					 ModelTimeExpCavesGDist> >;





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






template class FitOnlyRunner<System<ModelTimeExpCavesGDist,
				    ModelTimeExpCavesGDist>,
			     MyopicAgent<ModelTimeExpCavesGDist> >;

template class FitOnlyRunner<System<ModelTimeGDistTrendPow,
				    ModelTimeGDistTrendPow>,
			     MyopicAgent<ModelTimeGDistTrendPow> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDistTrendPowCon,
				    ModelTimeExpCavesGDistTrendPowCon>,
			     MyopicAgent<ModelTimeExpCavesGDistTrendPowCon> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDistTrendPowCon,
				    ModelTimeExpCavesEDist>,
			     MyopicAgent<ModelTimeExpCavesEDist> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDistTrendPowCon,
				    ModelTimeExpCavesGDist>,
			     MyopicAgent<ModelTimeExpCavesGDist> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDist,
				    ModelTimeExpCavesEDist>,
			     MyopicAgent<ModelTimeExpCavesEDist> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDist,
				    ModelRadius>,
			     MyopicAgent<ModelRadius> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDist,
				    ModelGDist>,
			     MyopicAgent<ModelGDist> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDist,
				    ModelCovar>,
			     MyopicAgent<ModelCovar> >;

template class FitOnlyRunner<System<ModelTimeExpCavesGDist,
				    ModelGDistKern>,
			     MyopicAgent<ModelGDistKern> >;



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
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    RankAgent<ToyFeatures0<ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelTimeExpCavesGDist>,
		      RankAgent<ToyFeatures0<ModelTimeExpCavesGDist>,
				ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    RankAgent<ToyFeatures0<ModelGDist>,
		      ModelGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDist>,
		      RankAgent<ToyFeatures0<ModelGDist>,
				ModelGDist>,
		      ModelGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    RankAgent<ToyFeatures4<ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelTimeExpCavesGDist>,
		      RankAgent<ToyFeatures4<ModelTimeExpCavesGDist>,
				ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesEDist>,
	    RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
		      ModelTimeExpCavesEDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelTimeExpCavesEDist>,
		      RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
				ModelTimeExpCavesEDist>,
		      ModelTimeExpCavesEDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelRadius>,
	    RankAgent<ToyFeatures4<ModelRadius>,
		      ModelRadius>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelRadius>,
		      RankAgent<ToyFeatures4<ModelRadius>,
				ModelRadius>,
		      ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    RankAgent<ToyFeatures4<ModelGDist>,
		      ModelGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDist>,
		      RankAgent<ToyFeatures4<ModelGDist>,
				ModelGDist>,
		      ModelGDist> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDistKern>,
	    RankAgent<ToyFeatures4<ModelGDistKern>,
		      ModelGDistKern>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDistKern>,
		      RankAgent<ToyFeatures4<ModelGDistKern>,
				ModelGDistKern>,
		      ModelGDistKern> >;




template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelTimeExpCavesGDist>,
		      RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
				ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeGDistTrendPow,
		   ModelTimeGDistTrendPow>,
	    RankAgent<ToyFeatures5<ModelTimeGDistTrendPow>,
		      ModelTimeGDistTrendPow>,
	    M1SpOptim<System<ModelTimeGDistTrendPow,
			     ModelTimeGDistTrendPow>,
		      RankAgent<ToyFeatures5<ModelTimeGDistTrendPow>,
				ModelTimeGDistTrendPow>,
		      ModelTimeGDistTrendPow> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    RankAgent<ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon>,
	    M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			     ModelTimeExpCavesGDistTrendPowCon>,
		      RankAgent<ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
				ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    RankAgent<ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon>,
	    M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			     ModelTimeExpCavesGDistTrendPowCon>,
		      RankAgent<ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>,
				ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    RankAgent<ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon>,
	    M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			     ModelTimeExpCavesGDistTrendPowCon>,
		      RankAgent<ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>,
				ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesEDist>,
	    RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
		      ModelTimeExpCavesEDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			     ModelTimeExpCavesEDist>,
		      RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
				ModelTimeExpCavesEDist>,
		      ModelTimeExpCavesEDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDist>,
	    RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			     ModelTimeExpCavesGDist>,
		      RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
				ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelRadius>,
	    RankAgent<ToyFeatures5<ModelRadius>,
		      ModelRadius>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelRadius>,
		      RankAgent<ToyFeatures5<ModelRadius>,
				ModelRadius>,
		      ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    RankAgent<ToyFeatures5<ModelGDist>,
		      ModelGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDist>,
		      RankAgent<ToyFeatures5<ModelGDist>,
				ModelGDist>,
		      ModelGDist> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelCovar>,
	    RankAgent<ToyFeatures5<ModelCovar>,
		      ModelCovar>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelCovar>,
		      RankAgent<ToyFeatures5<ModelCovar>,
				ModelCovar>,
		      ModelCovar> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDistKern>,
	    RankAgent<ToyFeatures5<ModelGDistKern>,
		      ModelGDistKern>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDistKern>,
		      RankAgent<ToyFeatures5<ModelGDistKern>,
				ModelGDistKern>,
		      ModelGDistKern> >;






template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    RankAgent<WnsFeatures1<ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelTimeExpCavesGDist>,
		      RankAgent<WnsFeatures1<ModelTimeExpCavesGDist>,
				ModelTimeExpCavesGDist>,
		      ModelTimeExpCavesGDist> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelRadius>,
	    RankAgent<WnsFeatures1<ModelRadius>,
		      ModelRadius>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelRadius>,
		      RankAgent<WnsFeatures1<ModelRadius>,
				ModelRadius>,
		      ModelRadius> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    RankAgent<WnsFeatures1<ModelGDist>,
		      ModelGDist>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDist>,
		      RankAgent<WnsFeatures1<ModelGDist>,
				ModelGDist>,
		      ModelGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDistKern>,
	    RankAgent<WnsFeatures1<ModelGDistKern>,
		      ModelGDistKern>,
	    M1SpOptim<System<ModelTimeExpCavesGDist,
			     ModelGDistKern>,
		      RankAgent<WnsFeatures1<ModelGDistKern>,
				ModelGDistKern>,
		      ModelGDistKern> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    OsspAgent<ModelTimeExpCavesGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelTimeExpCavesGDist>,
			OsspAgent<ModelTimeExpCavesGDist>,
			ToyFeatures0<ModelTimeExpCavesGDist>,
			ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    OsspAgent<ModelGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelGDist>,
			OsspAgent<ModelGDist>,
			ToyFeatures0<ModelGDist>,
			ModelGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    OsspAgent<ModelTimeExpCavesGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelTimeExpCavesGDist>,
			OsspAgent<ModelTimeExpCavesGDist>,
			ToyFeatures4<ModelTimeExpCavesGDist>,
			ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesEDist>,
	    OsspAgent<ModelTimeExpCavesEDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelTimeExpCavesEDist>,
			OsspAgent<ModelTimeExpCavesEDist>,
			ToyFeatures5<ModelTimeExpCavesEDist>,
			ModelTimeExpCavesEDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelRadius>,
	    OsspAgent<ModelRadius>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelRadius>,
			OsspAgent<ModelRadius>,
			ToyFeatures4<ModelRadius>,
			ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    OsspAgent<ModelGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelGDist>,
			OsspAgent<ModelGDist>,
			ToyFeatures4<ModelGDist>,
			ModelGDist> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    OsspAgent<ModelTimeExpCavesGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelTimeExpCavesGDist>,
			OsspAgent<ModelTimeExpCavesGDist>,
			ToyFeatures5<ModelTimeExpCavesGDist>,
			ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeGDistTrendPow,
		   ModelTimeGDistTrendPow>,
	    OsspAgent<ModelTimeGDistTrendPow>,
	    M1OsspOptim<System<ModelTimeGDistTrendPow,
			       ModelTimeGDistTrendPow>,
			OsspAgent<ModelTimeGDistTrendPow>,
			ToyFeatures5<ModelTimeGDistTrendPow>,
			ModelTimeGDistTrendPow> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
	    M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			       ModelTimeExpCavesGDistTrendPowCon>,
			OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
			ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
			ModelTimeExpCavesGDistTrendPowCon> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
	    M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			       ModelTimeExpCavesGDistTrendPowCon>,
			OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
			ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>,
			ModelTimeExpCavesGDistTrendPowCon> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDistTrendPowCon>,
	    OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
	    M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			       ModelTimeExpCavesGDistTrendPowCon>,
			OsspAgent<ModelTimeExpCavesGDistTrendPowCon>,
			ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>,
			ModelTimeExpCavesGDistTrendPowCon> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesEDist>,
	    OsspAgent<ModelTimeExpCavesEDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			       ModelTimeExpCavesEDist>,
			OsspAgent<ModelTimeExpCavesEDist>,
			ToyFeatures5<ModelTimeExpCavesEDist>,
			ModelTimeExpCavesEDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		   ModelTimeExpCavesGDist>,
	    OsspAgent<ModelTimeExpCavesGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			       ModelTimeExpCavesGDist>,
			OsspAgent<ModelTimeExpCavesGDist>,
			ToyFeatures5<ModelTimeExpCavesGDist>,
			ModelTimeExpCavesGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelRadius>,
	    OsspAgent<ModelRadius>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelRadius>,
			OsspAgent<ModelRadius>,
			ToyFeatures5<ModelRadius>,
			ModelRadius> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    OsspAgent<ModelGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelGDist>,
			OsspAgent<ModelGDist>,
			ToyFeatures5<ModelGDist>,
			ModelGDist> >;


template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelCovar>,
	    OsspAgent<ModelCovar>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelCovar>,
			OsspAgent<ModelCovar>,
			ToyFeatures5<ModelCovar>,
			ModelCovar> >;



template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelTimeExpCavesGDist>,
	    OsspAgent<ModelTimeExpCavesGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelTimeExpCavesGDist>,
			OsspAgent<ModelTimeExpCavesGDist>,
			WnsFeatures1<ModelTimeExpCavesGDist>,
			ModelTimeExpCavesGDist> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelRadius>,
	    OsspAgent<ModelRadius>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelRadius>,
			OsspAgent<ModelRadius>,
			WnsFeatures1<ModelRadius>,
			ModelRadius> >;

template class
OptimRunner<System<ModelTimeExpCavesGDist,
		   ModelGDist>,
	    OsspAgent<ModelGDist>,
	    M1OsspOptim<System<ModelTimeExpCavesGDist,
			       ModelGDist>,
			OsspAgent<ModelGDist>,
			WnsFeatures1<ModelGDist>,
			ModelGDist> >;



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
TuneRunner<System<ModelTimeExpCavesGDist,
		  ModelTimeExpCavesGDist>,
	   RankAgent<ToyFeatures0<ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist>,
	   M1SpOptim<System<ModelTimeExpCavesGDist,
			    ModelTimeExpCavesGDist>,
		     RankAgent<ToyFeatures0<ModelTimeExpCavesGDist>,
			       ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist> >;


template class
TuneRunner<System<ModelTimeExpCavesEDist,
		  ModelTimeExpCavesEDist>,
	   RankAgent<ToyFeatures0<ModelTimeExpCavesEDist>,
		     ModelTimeExpCavesEDist>,
	   M1SpOptim<System<ModelTimeExpCavesEDist,
			    ModelTimeExpCavesEDist>,
		     RankAgent<ToyFeatures0<ModelTimeExpCavesEDist>,
			       ModelTimeExpCavesEDist>,
		     ModelTimeExpCavesEDist> >;


template class
TuneRunner<System<ModelGDist,
		  ModelGDist>,
	   RankAgent<ToyFeatures0<ModelGDist>,
		     ModelGDist>,
	   M1SpOptim<System<ModelGDist,
			    ModelGDist>,
		     RankAgent<ToyFeatures0<ModelGDist>,
			       ModelGDist>,
		     ModelGDist> >;


template class
TuneRunner<System<ModelTimeExpCavesGDist,
		  ModelTimeExpCavesGDist>,
	   RankAgent<ToyFeatures4<ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist>,
	   M1SpOptim<System<ModelTimeExpCavesGDist,
			    ModelTimeExpCavesGDist>,
		     RankAgent<ToyFeatures4<ModelTimeExpCavesGDist>,
			       ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist> >;


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
TuneRunner<System<ModelGDist,
		  ModelGDist>,
	   RankAgent<ToyFeatures4<ModelGDist>,
		     ModelGDist>,
	   M1SpOptim<System<ModelGDist,
			    ModelGDist>,
		     RankAgent<ToyFeatures4<ModelGDist>,
			       ModelGDist>,
		     ModelGDist> >;


template class
TuneRunner<System<ModelGDistKern,
		  ModelGDistKern>,
	   RankAgent<ToyFeatures4<ModelGDistKern>,
		     ModelGDistKern>,
	   M1SpOptim<System<ModelGDistKern,
			    ModelGDistKern>,
		     RankAgent<ToyFeatures4<ModelGDistKern>,
			       ModelGDistKern>,
		     ModelGDistKern> >;


template class
TuneRunner<System<ModelTimeExpCavesGDist,
		  ModelTimeExpCavesGDist>,
	   RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist>,
	   M1SpOptim<System<ModelTimeExpCavesGDist,
			    ModelTimeExpCavesGDist>,
		     RankAgent<ToyFeatures5<ModelTimeExpCavesGDist>,
			       ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist> >;


template class
TuneRunner<System<ModelTimeGDistTrendPow,
		  ModelTimeGDistTrendPow>,
	   RankAgent<ToyFeatures5<ModelTimeGDistTrendPow>,
		     ModelTimeGDistTrendPow>,
	   M1SpOptim<System<ModelTimeGDistTrendPow,
			    ModelTimeGDistTrendPow>,
		     RankAgent<ToyFeatures5<ModelTimeGDistTrendPow>,
			       ModelTimeGDistTrendPow>,
		     ModelTimeGDistTrendPow> >;


template class
TuneRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		  ModelTimeExpCavesGDistTrendPowCon>,
	   RankAgent<ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
		     ModelTimeExpCavesGDistTrendPowCon>,
	   M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			    ModelTimeExpCavesGDistTrendPowCon>,
		     RankAgent<ToyFeatures5<ModelTimeExpCavesGDistTrendPowCon>,
			       ModelTimeExpCavesGDistTrendPowCon>,
		     ModelTimeExpCavesGDistTrendPowCon> >;


template class
TuneRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		  ModelTimeExpCavesGDistTrendPowCon>,
	   RankAgent<ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>,
		     ModelTimeExpCavesGDistTrendPowCon>,
	   M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			    ModelTimeExpCavesGDistTrendPowCon>,
		     RankAgent<ToyFeatures6<ModelTimeExpCavesGDistTrendPowCon>,
			       ModelTimeExpCavesGDistTrendPowCon>,
		     ModelTimeExpCavesGDistTrendPowCon> >;


template class
TuneRunner<System<ModelTimeExpCavesGDistTrendPowCon,
		  ModelTimeExpCavesGDistTrendPowCon>,
	   RankAgent<ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>,
		     ModelTimeExpCavesGDistTrendPowCon>,
	   M1SpOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			    ModelTimeExpCavesGDistTrendPowCon>,
		     RankAgent<ToyFeatures7<ModelTimeExpCavesGDistTrendPowCon>,
			       ModelTimeExpCavesGDistTrendPowCon>,
		     ModelTimeExpCavesGDistTrendPowCon> >;


template class
TuneRunner<System<ModelTimeExpCavesEDist,
		  ModelTimeExpCavesEDist>,
	   RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
		     ModelTimeExpCavesEDist>,
	   M1SpOptim<System<ModelTimeExpCavesEDist,
			    ModelTimeExpCavesEDist>,
		     RankAgent<ToyFeatures5<ModelTimeExpCavesEDist>,
			       ModelTimeExpCavesEDist>,
		     ModelTimeExpCavesEDist> >;


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
TuneRunner<System<ModelGDist,
		  ModelGDist>,
	   RankAgent<ToyFeatures5<ModelGDist>,
		     ModelGDist>,
	   M1SpOptim<System<ModelGDist,
			    ModelGDist>,
		     RankAgent<ToyFeatures5<ModelGDist>,
			       ModelGDist>,
		     ModelGDist> >;


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
TuneRunner<System<ModelGDistKern,
		  ModelGDistKern>,
	   RankAgent<ToyFeatures5<ModelGDistKern>,
		     ModelGDistKern>,
	   M1SpOptim<System<ModelGDistKern,
			    ModelGDistKern>,
		     RankAgent<ToyFeatures5<ModelGDistKern>,
			       ModelGDistKern>,
		     ModelGDistKern> >;



template class
TuneRunner<System<ModelTimeExpCavesGDist,
		  ModelTimeExpCavesGDist>,
	   RankAgent<WnsFeatures1<ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist>,
	   M1SpOptim<System<ModelTimeExpCavesGDist,
			    ModelTimeExpCavesGDist>,
		     RankAgent<WnsFeatures1<ModelTimeExpCavesGDist>,
			       ModelTimeExpCavesGDist>,
		     ModelTimeExpCavesGDist> >;



template class
TuneRunner<System<ModelTimeExpCavesEDist,
		  ModelTimeExpCavesEDist>,
	   RankAgent<WnsFeatures1<ModelTimeExpCavesEDist>,
		     ModelTimeExpCavesEDist>,
	   M1SpOptim<System<ModelTimeExpCavesEDist,
			    ModelTimeExpCavesEDist>,
		     RankAgent<WnsFeatures1<ModelTimeExpCavesEDist>,
			       ModelTimeExpCavesEDist>,
		     ModelTimeExpCavesEDist> >;


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
TuneRunner<System<ModelGDist,
		  ModelGDist>,
	   RankAgent<WnsFeatures1<ModelGDist>,
		     ModelGDist>,
	   M1SpOptim<System<ModelGDist,
			    ModelGDist>,
		     RankAgent<WnsFeatures1<ModelGDist>,
			       ModelGDist>,
		     ModelGDist> >;




template class
TuneRunner<System<ModelGDistKern,
		  ModelGDistKern>,
	   RankAgent<WnsFeatures1<ModelGDistKern>,
		     ModelGDistKern>,
	   M1SpOptim<System<ModelGDistKern,
			    ModelGDistKern>,
		     RankAgent<WnsFeatures1<ModelGDistKern>,
			       ModelGDistKern>,
		     ModelGDistKern> >;


template class
TuneRunner<System<ModelGravityGDist,
		  ModelGravityGDist>,
	   RankAgent<WnsFeatures1<ModelGravityGDist>,
		     ModelGravityGDist>,
	   M1SpOptim<System<ModelGravityGDist,
			    ModelGravityGDist>,
		     RankAgent<WnsFeatures1<ModelGravityGDist>,
			       ModelGravityGDist>,
		     ModelGravityGDist> >;


template class
TuneRunner<System<ModelGravityEDist,
		  ModelGravityEDist>,
	   RankAgent<WnsFeatures1<ModelGravityEDist>,
		     ModelGravityEDist>,
	   M1SpOptim<System<ModelGravityEDist,
			    ModelGravityEDist>,
		     RankAgent<WnsFeatures1<ModelGravityEDist>,
			       ModelGravityEDist>,
		     ModelGravityEDist> >;






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
