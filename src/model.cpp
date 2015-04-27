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
  for(i = 0; i < numPars; ++i){
    pars[i]->setFill(probs,sD,tD,fD,dD);
  }
  set = 1;
  ready = 0;
}


void ModelBase::modFill(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  if(set == 1){
    int i,numPars = pars.size();
    for(i = 0; i < numPars; ++i){
      pars[i]->modFill(probs,sD,tD,fD,dD);
    }
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


std::vector<double> ModelBase::partial(const int notNode,
				       const int infNode,
				       const SimData & sD,
				       const TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD){
  std::vector<double> p,pi;
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    pi = pars[i]->partial(notNode,infNode,sD,tD,fD,dD);
    p.insert(p.end(),pi.begin(),pi.end());
  }
  return p;
}


std::vector<double> ModelBase::partial2(const int notNode,
					const int infNode,
					const SimData & sD,
					const TrtData & tD,
					const FixedData & fD,
					const DynamicData & dD){
  std::vector<double> p,pi;
  std::vector<int> parsLen;
  int i,numPars = pars.size(),totLen = 0;
  for(i = 0; i < numPars; ++i){
    parsLen.push_back(totLen);
    totLen += pars[i]->size();
  }
  parsLen.push_back(totLen);


  p.resize(totLen*totLen);
  std::fill(p.begin(),p.end(),0);

  int j,k,n,N,D;
  for(i = 0; i < numPars; ++i){
    pi = pars[i]->partial2(notNode,infNode,sD,tD,fD,dD);

    n = parsLen[i];
    N = parsLen[i+1];
    D = N - n;
    for(j = n; j < N; ++j){
      for(k = j; k < N; ++k){
	p[j*totLen + k] = pi[(j-n)*D + (k-n)];
	p[k*totLen + j] = pi[(k-n)*D + (j-n)];
      }
    }
    
  }
  
  return p;
}


void ModelBase::setFisher(const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
  int i,parsSize = pars.size(),totPar = 0;
  std::vector<int> parsLen;
  for(i = 0; i < parsSize; ++i){
    parsLen.push_back(totPar);
    totPar += pars[i]->size();
  }
  parsLen.push_back(totPar);

  fisher.resize(totPar*totPar);
  std::fill(fisher.begin(),fisher.end(),0);
  
  
  
  std::vector<std::vector<int> > hist = sD.history;
  hist.push_back(sD.status);
  std::vector<DataBundle> db = historyToData(hist);

  int t,iN,nN,pi,pj;
  for(t = 0; t < sD.time; ++t){
    SimData sDi= std::get<0>(db[t]);
    TrtData tDi= std::get<1>(db[t]);
    DynamicData dDi = std::get<2>(db[t]);

    setFill(sDi,tDi,fD,dDi);
    setQuick(sDi,tDi,fD,dDi);
    infProbs(sDi,tDi,fD,dDi);

    for(nN = 0; nN < sD.numNotInfec; ++nN){
      std::vector<double> dbl(totPar*totPar,0);
      std::vector<double> sqr(totPar,0);
      
      for(iN = 0; iN < sD.numInfected; ++iN){
	int ind = nN*sD.numInfected + iN;
	std::vector<double> p = partial(sD.notInfec[nN],
					sD.infected[iN],
					sDi,tDi,fD,dDi);
	std::vector<double> p2 = partial2(sD.notInfec[nN],
					  sD.infected[iN],
					  sDi,tDi,fD,dDi);

	for(pi = 0; pi < totPar; ++pi){
	  sqr[pi] += quick[ind]*p[pi];
	  
	  int piInd = pi*totPar;
	  for(pj = pi; pj < totPar; ++pj){
	    dbl[piInd + pj] += quick[ind]*p2[piInd + pj];
	    dbl[piInd + pj] += quick[ind]*p[pj]*p[pi];

	    dbl[pj*totPar + pi] += quick[ind]*p2[pj*totPar + pi];
	    dbl[pj*totPar + pi] += quick[ind]*p[pi]*p[pj];
	  }
	}
      }

      int next = (hist[t+1][sD.notInfec[nN]] < 2 ? 0 : 1);

      for(pi = 0; pi < totPar; ++pi){
	for(pj = pi; pj < totPar; ++pj){
	  double prob = expitInfProbs[nN];
	  
	  fisher[pi*totPar + pj] +=
	    (double(next)/prob - 1.0)*dbl[pi*totPar + pj];
	  
	  fisher[pj*totPar + pi] +=
	    (double(next)/prob - 1.0)*dbl[pj*totPar + pi];
	  

	  if(next == 1){
	    fisher[pi*totPar + pj] -=
	      (1-prob)*sqr[pi]*sqr[pj]/(prob*prob);
	    
	    fisher[pj*totPar + pi] -=
	      (1-prob)*sqr[pj]*sqr[pi]/(prob*prob);
	  }
	}
      }
    }
  }
}
