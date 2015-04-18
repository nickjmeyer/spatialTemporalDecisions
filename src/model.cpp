#include "model.hpp"


ModelBase::ModelBase(const std::vector<ParamBase *> & newPars,
		     const FixedData & fD){
  pars = newPars;
  std::for_each(pars.begin(),pars.end(),
		[&fD](ParamBase * p){
		  p->init(fD);
		});
}


ModelBase::~ModelBase(){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    delete pars.at(i);
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
  expitInfProbs.resize(sD.numNotInfec);
  int i,j,node0;
  double prob;
  for(i=0; i<sD.numNotInfec; i++){
    node0 = sD.notInfec.at(i) * fD.numNodes;
    prob=1.0;
    for(j=0; j<sD.numInfected; j++)
      prob *= 1.0/(1.0+std::exp(probs.at(node0 + sD.infected.at(j))));
    expitInfProbs.at(i) = 1.0-prob;
  }
}


std::vector<double> ModelBase::infProbs(){
  return expitInfProbs;
}


void ModelBase::revProbs(const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD){
  expitRevProbs.clear();
  expitRevProbs.reserve(sD.numInfected);
  int i,j,node0;
  double prob;
  for(i=0; i<sD.numInfected; i++){
    node0 = sD.infected.at(i);
    prob=1.0;
    for(j=0; j<sD.numNotInfec; j++)
      prob *= 1.0/(1.0+std::exp(probs.at(node0 +
					 sD.notInfec.at(j)*fD.numNodes)));
    expitRevProbs.push_back(1.0-prob);
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
    pars.at(i)->setFill(probs,sD,tD,fD,dD);
  set = 1;
}


void ModelBase::modFill(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  if(set == 1){
    int i,numPars = pars.size();
    for(i = 0; i < numPars; ++i)
      pars.at(i)->modFill(probs,sD,tD,fD,dD);
  }
  else if(set == 0){
    setFill(sD,tD,fD,dD);
  }
}



double ModelBase::oneOnOne(const int notNode,
			   const int infNode,
			   const int numNodes) const {
  return probs.at(notNode * numNodes + infNode);
}


std::vector<double> ModelBase::getPar() const{
  std::vector<double> all,add;
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    add = pars.at(i)->getPar();
    all.insert(all.end(),add.begin(),add.end());
  }
  return all;
}


std::vector<double>::const_iterator
ModelBase::putPar(std::vector<double>::const_iterator it){
  set = 0;
  
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    it = pars.at(i)->putPar(it);
  }
  return it;
}
