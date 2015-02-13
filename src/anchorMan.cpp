#include "anchorMan.hpp"


std::vector<double> AnchorManTunePar::getPar() const{
  return std::vector<double>();
}


void AnchorManTunePar::putPar(const std::vector<double> & par){
  // do nothing;
}


AnchorManTunePar::AnchorManTunePar(){
  numSamples = 100;
  cutoff = .95;
  freq = 6;
}

template class AnchorMan<System<GravityModel,GravityParam,
				GravityModel,GravityParam>,
			 RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
				      GravityModel,GravityParam>,
			 FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
				     GravityModel,GravityParam>,
			 GravityModel,GravityParam>;

template class AnchorMan<System<RangeModel,RangeParam,
				RangeModel,RangeParam>,
			 RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				      RangeModel,RangeParam>,
			 FeaturesInt<ToyFeatures2<RangeModel,RangeParam>,
				     RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class AnchorMan<System<GravityModel,GravityParam,
				RangeModel,RangeParam>,
			 RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				      RangeModel,RangeParam>,
			 FeaturesInt<ToyFeatures2<RangeModel,RangeParam>,
				     RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class AnchorMan<System<GravityModel,GravityParam,
				CaveModel,CaveParam>,
			 RankToyAgent<ToyFeatures2<CaveModel,CaveParam>,
				      CaveModel,CaveParam>,
			 FeaturesInt<ToyFeatures2<CaveModel,CaveParam>,
				     CaveModel,CaveParam>,
			 CaveModel,CaveParam>;


template <class S, class A, class F,
	  class M, class MP>
AnchorMan<S,A,F,M,MP>::AnchorMan(){
  // set switched time to max
  name = "AM";
  switched = std::numeric_limits<int>::max();
}


template <class S, class A, class F,
	  class M, class MP>
void AnchorMan<S,A,F,M,MP> 
::optim(const S & system, A & agent){

  if(switched == std::numeric_limits<int>::max())
    switched = system.sD.time;

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);
  
  std::vector<double> w = agent.tp.getPar();

  std::string msg;
  msg = "vals: (" + njm::toString(system.sD.time,"",0,0)
    + ", " + njm::toString(switched,"",0,0) + ")\n";

  if(system.sD.time == system.fD.trtStart ||
     (system.sD.time - system.fD.trtStart) % tp.freq != 0){ // do normal thing
    msg += " -- decision: normal";
    if(switched >= system.sD.time){
      msg += " -- optim: m1";
      m1Opt.optim(s,agent);
      switched+=system.fD.period;
    }
    else{
      msg += " -- optim: m2";
      m2Opt.optim(s,agent);
    }
  }
  else{ // check switch if not already switched
    msg += " -- decision: Anchor Man";
    if(switched >= system.sD.time){ // haven't made switch yet
      msg += " -- switched: no";
  
      m1Opt.optim(s,agent);
      m1W = agent.tp.getPar();

      agent.tp.putPar(w);
      m2Opt.optim(s,agent);
      m2W = agent.tp.getPar();

      // get ready for generating null distribution
      s = System<M,MP,M,MP>(system.sD_r,system.tD_r,system.fD,system.dD_r,
			    system.modelEst,system.modelEst,
			    system.paramEst,system.paramEst_r);
    
      if(toSwitch(s,agent,system.sD.time)){
	msg += " -- toswitch: yes";
	agent.tp.putPar(m2W);
	switched = system.sD.time;
      }
      else{
	msg += " -- toswitch: no";
	agent.tp.putPar(m1W);
	switched+=system.fD.period;
      }
    }
    else{ // made switch, use m2
      msg += " -- switched: yes";
      m2Opt.optim(s,agent);
    }
  }
  njm::message(msg);
}


template <class S, class A, class F,
	  class M, class MP>
int AnchorMan<S,A,F,M,MP>
::toSwitch(System<M,MP,M,MP> & system, A & agent, const int T){
  double testStat = njm::l2norm(m1W,m2W);

  std::vector<double> samples;
  samples.reserve(tp.numSamples);
  int i, t;
  for(i = 0; i < tp.numSamples; i++){
    system.reset();
    for(t = 0; t < T; t++){
      if(omp_get_thread_num() == 0)
	std::cout << "\r(" << i << ", " << T << ", "
		  << system.sD.time << ", " << t << ")        "
		  << std::flush;
      if(t >= system.fD.trtStart &&
	 (((t - system.fD.trtStart) % system.fD.period) == 0)){
	system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
			    system.paramEst);
	m1Opt.optim(system,agent);
      }
      
      if(t >= system.fD.trtStart)
	agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
		       system.modelEst,system.paramEst);
      system.updateStatus();
      system.nextPoint();
    }

    if(omp_get_thread_num() == 0)
      std::cout << "\r(" << i << ", " << T << ", "
		<< system.sD.time << ", " << t << ")        "
		<< std::flush;

    m2Opt.optim(system,agent);

    samples.push_back(njm::l2norm(m1W,agent.tp.getPar()));
  }

  std::sort(samples.begin(),samples.end());

  double mean,sd,mn,mx;
  mean = 0;
  sd = 0;
  mn = std::numeric_limits<double>::max();
  mx = std::numeric_limits<double>::lowest();
  for(i = 0; i < tp.numSamples; i++){
    if(samples.at(i) < mn)
      mn = samples.at(i);
    if(samples.at(i) > mx)
      mx = samples.at(i);
    mean+=samples.at(i);
    sd+=samples.at(i)*samples.at(i);
  }

  mean/=double(tp.numSamples);
  sd /=double(tp.numSamples);
  sd -= mean*mean;

  int ind = int(std::ceil(tp.cutoff*tp.numSamples + 0.5)) - 1;

  if(omp_get_thread_num() == 0){
    std::cout << std::endl
	      << "test: " << testStat << " > " << samples.at(ind)
	      << "[" << mean << ", " << sd << "]"
	      << std::endl;
  }
  

  return (testStat > samples.at(ind) ? 1 : 0);
}

