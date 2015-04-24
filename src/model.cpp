#include "model.hpp"


ModelBase::ModelBase(const std::vector<ParamBase *> & newPars,
		     const FixedData & fD){
  set = 0;
  ready = 0;
  pars = newPars;
  std::for_each(pars.begin(),pars.end(),
		[&fD](ParamBase * p){
		  p->init(fD);
		});
}


ModelBase::~ModelBase(){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    delete pars[i];
  }
}


void ModelBase::setType(const Estimation & est){
  fitType = est;
}


Estimation ModelBase::getType() const{
  return fitType;
}


Estimation & ModelBase::getType() {
  return fitType;
}


void ModelBase::infProbs(const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD){
  if(ready == 1){
    expitInfProbs.resize(sD.numNotInfec);
    int i,j,k;
    double prob;
    for(i = 0, k = 0; i < sD.numNotInfec; ++i){
      prob=1.0;
      for(j = 0; j < sD.numInfected; ++j,++k)
	prob *= quick[k];
      expitInfProbs[i] = 1.0-prob;
    }
  }
  else if(ready == 0){
    expitInfProbs.resize(sD.numNotInfec);
    int i,j,k;
    double prob;
    for(i = 0; i < sD.numNotInfec; ++i){
      k = sD.notInfec[i] * fD.numNodes;
      prob=1.0;
      for(j = 0; j < sD.numInfected; ++j)
	prob *= 1.0 / (1.0 + std::exp(probs[k + sD.infected[j]]));
      expitInfProbs[i] = 1.0-prob;
    }
  }    
}


std::vector<double> ModelBase::infProbs(){
  return expitInfProbs;
}


void ModelBase::revProbs(const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD){
  if(ready == 1){
    expitRevProbs.resize(sD.numInfected);
    int i,j,k;
    double prob;
    for(i = 0,k = 0; i < sD.numInfected; ++i){
      prob=1.0;
      for(j = 0; j < sD.numNotInfec; ++j,++k)
	prob *= quick[k];
      expitRevProbs[i] = 1.0-prob;
    }
  }
  else if(ready == 0){
    expitRevProbs.resize(sD.numInfected);
    int i,j,k;
    double prob;
    for(i = 0; i < sD.numInfected; ++i){
      k = sD.infected[i];
      prob=1.0;
      for(j = 0; j < sD.numNotInfec; ++j)
	prob *= 1.0 / (1.0 + std::exp(probs[k + sD.notInfec[j]*fD.numNodes]));
      expitRevProbs[i] = 1.0-prob;
    }
  }
  else{
    throw(1);
  }
}


std::vector<double> ModelBase::revProbs(){
  return expitRevProbs;
}



void ModelBase::setFill(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  int i,numPars = pars.size();
  probs = std::vector<double>(fD.numNodes*fD.numNodes,0.0);
  for(i = 0; i < numPars; ++i)
    pars[i]->setFill(probs,sD,tD,fD,dD);
  set = 1;
  ready = 0;
}


void ModelBase::modFill(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  if(set == 1){
    int i,numPars = pars.size();
    for(i = 0; i < numPars; ++i)
      pars[i]->modFill(probs,sD,tD,fD,dD);
  }
  else if(set == 0){
    setFill(sD,tD,fD,dD);
  }
  ready = 0;
}


void ModelBase::setQuick(const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD){
  int i,j,k,pK;
  quick.resize(sD.numNotInfec * sD.numInfected);
  for(i = 0,k = 0; i < sD.numNotInfec; ++i){
    pK = sD.notInfec[i]*fD.numNodes;
    for(j = 0; j < sD.numInfected; ++j,++k){
      quick[k] = 1.0/(1.0 + std::exp(probs[pK + sD.infected[j]]));
    }
  }
}


std::vector<double> & ModelBase::getQuick() {
  return quick;
}



double ModelBase::oneOnOne(const int notNode,
			   const int infNode,
			   const int numNodes) const {
  return probs[notNode * numNodes + infNode];
}



std::vector<double> ModelBase::getPar() const{
  std::vector<double> all,add;
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    add = pars[i]->getPar();
    all.insert(all.end(),add.begin(),add.end());
  }
  return all;
}


std::vector<double>::const_iterator
ModelBase::putPar(std::vector<double>::const_iterator it){
  set = 0;
  
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    it = pars[i]->putPar(it);
  }
  return it;
}
