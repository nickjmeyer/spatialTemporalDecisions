#include "modelIntercept.hpp"


static std::vector<ParamBase *> genPars(){
  std::vector<ParamBase *> pars;
  pars.push_back(new ParamIntercept);
  pars.push_back(new ParamTrt);
  return pars;
}

ModelIntercept::ModelIntercept(const FixedData & fD)
  : ModelBase("Intercept",genPars(),fD){
}


ModelIntercept::ModelIntercept(const ModelIntercept & m){
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
  fixSample = m.fixSample;
  setEdgeToEdge(m.getEdgeToEdge());
}


ModelIntercept & ModelIntercept::operator=(const ModelIntercept & m){
  if(this != & m){
    this->ModelIntercept::~ModelIntercept();
    new (this) ModelIntercept(m);
  }
  return *this;
}
