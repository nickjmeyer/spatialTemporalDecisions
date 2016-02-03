#include "model2GravityGDist.hpp"

static std::vector<ParamBase *> genPars(){
  std::vector<ParamBase *> pars;
  pars.push_back(new ParamIntercept);
  pars.push_back(new ParamBeta2);
  pars.push_back(new ParamGravityGDist);
  pars.push_back(new ParamTrt);
  return pars;
}

Model2GravityGDist::Model2GravityGDist(const FixedData & fD)
  : ModelBase("2GravityGDist",genPars(),fD){
}


Model2GravityGDist::Model2GravityGDist(const Model2GravityGDist & m){
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
  meanHit = m.meanHit;
  varHit = m.varHit;
  ready = m.ready;
  numInfected = m.numInfected;
  numNotInfec = m.numNotInfec;
  fitType = m.fitType;
  mcmc = m.mcmc;
  fixSample = m.fixSample;
}


Model2GravityGDist &
Model2GravityGDist::operator=(const Model2GravityGDist & m){
  if(this != & m){
    this->Model2GravityGDist::~Model2GravityGDist();
    new (this) Model2GravityGDist(m);
  }
  return *this;
}
