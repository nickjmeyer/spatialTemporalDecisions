#include "modelGDist.hpp"


static std::vector<ParamBase *> genPars(){
  std::vector<ParamBase *> pars;
  pars.push_back(new ParamIntercept);
  pars.push_back(new ParamGDist);
  pars.push_back(new ParamTrt);
  return pars;
}

ModelGDist::ModelGDist(const FixedData & fD)
  : ModelBase("GDist",genPars(),fD){
}


ModelGDist::ModelGDist(const ModelGDist & m){
  int i, parsSize = m.pars.size();
  pars.clear();
  for(i = 0; i < parsSize; ++i)
    pars.push_back(m.pars.at(i)->clone());

  name = m.name;
  numPars = m.numPars;
  set = m.set;
  probs = m.probs;
  expitInfProbs = m.expitInfProbs;
  expitRevProbs = m.expitRevProbs;
  quick = m.quick;
  fisher = m.fisher;
  meanHit = m.meanHit;
  varHit = m.varHit;
  ready = m.ready;
  numInfected = m.numInfected;
  numNotInfec = m.numNotInfec;
  fitType = m.fitType;
  mcmc = m.mcmc;
  fixSample = m.fixSample;
}


ModelGDist & ModelGDist::operator=(const ModelGDist & m){
  if(this != & m){
    this->ModelGDist::~ModelGDist();
    new (this) ModelGDist(m);
  }
  return *this;
}
