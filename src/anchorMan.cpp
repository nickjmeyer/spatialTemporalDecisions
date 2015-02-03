#include "anchorMan.hpp"


std::vector<double> AnchorManTunePar::getPar() const{
  return std::vector<double>();
}


void AnchorManTunePar(const std::vector<double> & par){
  // do nothing;
}


AnchorManTunePar::AnchorManTunePar(){
  numSamples = 100;
  cutoff = .95;
}


template <class System, class Agent, class Features,
	  class Model, class ModelParam>
void AnchorMan<System,Agent,Features,Model,ModelParam>
::optim(System system, Agent & agent){
  if(toSwitch(system,agent))
    m1Opt.optim(system,agent);
  else
    m2Opt.optim(system,agent);
}


template <class System, class Agent, class Features,
	  class Model, class ModelParam>
int AnchorMan<System,Agent,Features,Model,ModelParam>
::toSwitch(System system, Agent & agent){
  std::vector<double> samples;
  samples.reserve(tp.numSamples);
  int i;
  for(i = 0; i < tp.numSamples; i++)
    samples.push_back(sampleNull(system,agent,system.sD.time));

  int ind = std::ceil(tp.cutoff*tp.numSamples + 0.5) - 1;

  double testStat = getTestStat(system,agent);

  return (testStat > samples.at(ind) ? 1 : 0);
}


template <class System, class Agent, class Features,
	  class Model, class ModelParam>
double AnchorMan<System,Agent,Features,Model,ModelParam>
::sampleNull(System system, Agent & agent,
	     const int numYears){

  system.reset();
  int t;
  for(t = 0; t < numYears; t++){
    if(t >= system.fD.trtStart &&
       (((t - system.fD.trtStart) % system.fD.period) == 0)){
      system.modelEst.fit(system.sD.system.tD,system.fD,system.dD,
			  system.paramEst);
      m1Opt.optim(system,agent);
    }

    if(t >= system.fD.trtStart)
      agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		     system.modelEstl,system.paramEst);
    system.updateStatus();
    system.nextPoint();
  }

  return getTestStat(system,agent);
  
}


template <class System, class Agent, class Features,
	  class Model, class ModelParam>
double AnchorMan<System,Agent,Features,Model,ModelParam>
::getTestStat(System system,Agent & agent){
  m2Opt.qEval.preCompData(system.sD,system.fD);
  m2Opt.qEval.bellResFixData(system.sD,system.tD,system.fD,system.dD,
			     system.modelEst,system.paramEst);
  m2Opt.qEval.bellResPolData(system.sD.time,system.fD,
			     system.modelEst,system.paramEst,agent);
  m2Opt.qEval.solve();
  return m2Opt.qFn(system.sD,system.tD,system.fD,system.dD,
		   system.modelEst,system.paramEst,agent);
}
