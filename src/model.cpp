#include "model.hpp"


ModelBase::ModelBase(const std::string & str,
		     const std::vector<ParamBase *> & newPars,
		     const FixedData & fD){
  name = str;
  
  set = 0;
  ready = 0;
  pars = newPars;
  std::for_each(pars.begin(),pars.end(),
		[&fD](ParamBase * p){
		  p->init(fD);
		});
  numPars = 0;
  std::for_each(pars.begin(),pars.end(),
		[this](ParamBase * p){
		  numPars += p->size();
		});
}


ModelBase::~ModelBase(){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    delete pars.at(i);
  }
}


void ModelBase::read(){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    pars[i]->read(name);
  }
}


void ModelBase::save() const{
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    pars[i]->save(name);
  }
}


void ModelBase::linScale(const double & scale){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i)
    pars.at(i)->linScale(scale);
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


std::vector<double>
ModelBase::getPar(const std::vector<std::string> & name) const{
  std::vector<double> res,add;
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i){
    add = pars[i]->getPar(name);
    res.insert(res.end(),add.begin(),add.end());
  }
  return res;
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


void ModelBase::setPar(const std::string & name,
		       const double & val){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i)
    pars[i]->setPar(name,val);
}


void ModelBase::setPar(const std::vector<std::string> & name,
		       const double & val){
  int i,numPars = pars.size();
  for(i = 0; i < numPars; ++i)
    pars[i]->setPar(name,val);
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
  fisher.resize(numPars * numPars);
  std::fill(fisher.begin(),fisher.end(),0);
  
  
  
  std::vector<std::vector<int> > hist = sD.history;
  hist.push_back(sD.status);
  std::vector<DataBundle> db = historyToData(hist);

  int t,iN,nN;
  unsigned int pi,pj;
  for(t = 0; t < sD.time; ++t){
    SimData sDi= std::get<0>(db[t]);
    TrtData tDi= std::get<1>(db[t]);
    DynamicData dDi = std::get<2>(db[t]);

    setFill(sDi,tDi,fD,dDi);
    setQuick(sDi,tDi,fD,dDi);
    infProbs(sDi,tDi,fD,dDi);

    for(nN = 0; nN < sDi.numNotInfec; ++nN){
      std::vector<double> dbl(numPars*numPars,0);
      std::vector<double> sqr(numPars,0);
      
      for(iN = 0; iN < sDi.numInfected; ++iN){
	int ind = nN*sDi.numInfected + iN;
	std::vector<double> p = partial(sDi.notInfec[nN],
					sDi.infected[iN],
					sDi,tDi,fD,dDi);
	std::vector<double> p2 = partial2(sDi.notInfec[nN],
					  sDi.infected[iN],
					  sDi,tDi,fD,dDi);
	
	double quickInd = std::max(1e-10,quick[ind]);
	for(pi = 0; pi < numPars; ++pi){
	  sqr[pi] += quickInd*p[pi];
	  
	  int piInd = pi*numPars;
	  for(pj = pi; pj < numPars; ++pj){
	    dbl[piInd + pj] += quickInd*p2[piInd + pj];
	    dbl[piInd + pj] += quickInd*p[pj]*p[pi];
	    if(pj != pi){
	      dbl[pj*numPars + pi] += quickInd*p2[pj*numPars + pi];
	      dbl[pj*numPars + pi] += quickInd*p[pi]*p[pj];
	    }
	  }
	}
      }

      int next = (hist[t+1][sDi.notInfec[nN]] < 2 ? 0 : 1);

      for(pi = 0; pi < numPars; ++pi){
	for(pj = pi; pj < numPars; ++pj){
	  double prob = std::max(1e-10,expitInfProbs[nN]);

	  fisher[pi*numPars + pj] +=
	    (double(next)/prob - 1.0)*dbl[pi*numPars + pj];

	  if(pj != pi){
	    fisher[pj*numPars + pi] +=
	      (double(next)/prob - 1.0)*dbl[pj*numPars + pi];
	  }
	  

	  if(next == 1){
	    fisher[pi*numPars + pj] -=
	      (1-prob)*sqr[pi]*sqr[pj]/(prob*prob);

	    if(pj != pi){
	      fisher[pj*numPars + pi] -=
		(1-prob)*sqr[pj]*sqr[pi]/(prob*prob);
	    }
	  }
	}
      }
    }
  }

  Eigen::Map<Eigen::MatrixXd> I(fisher.data(),numPars,numPars);
  if(sD.time <= fD.trtStart){
    I.row(numPars-2).setZero();
    I.row(numPars-1).setZero();
    I.col(numPars-2).setZero();
    I.col(numPars-1).setZero();

    I(numPars-2,numPars-2) = -2.0;
    I(numPars-1,numPars-1) = -2.0;
  }
  // std::cout << "I: " << std::endl
  // 	    << I << std::endl
  // 	    << std::endl;
  Eigen::LDLT<Eigen::MatrixXd> ldlt(-I);

  if(ldlt.info() != Eigen::Success){
    static int badCnt = 0;
    std::cout << "Error: ModelBase::getFisher(): eigen solver failed"
	      << std::endl;
    njm::toFile(fisher,njm::sett.datExt("badFisher" +
					njm::toString(badCnt++,"",0,0)));
    throw(1);
  }
  else{
    Eigen::Transpositions<Eigen::Dynamic,Eigen::Dynamic> P;
    P = ldlt.transpositionsP();
    Eigen::MatrixXd L = ldlt.matrixL();
    // std::cout << "L: " << std::endl
    // 	      << L << std::endl
    // 	      << std::endl;
    L.transposeInPlace();    
    Eigen::VectorXd D = ldlt.vectorD();
    // std::cout << "D: " << std::endl
    // 	      << D << std::endl
    // 	      << std::endl;
    varHit = P.transpose() *
      L.inverse() *
      D.cwiseInverse().cwiseAbs().cwiseSqrt().asDiagonal();
    std::vector<double> all = getPar();
    meanHit = Eigen::Map<Eigen::VectorXd>(all.data(),numPars);
  }
}


bool ModelBase::sample(){
  if(fitType == MLES){
    std::vector<double> rand;
    unsigned int i;
    for(i = 0; i < numPars; ++i){
      rand.push_back(njm::rnorm01());
    }
    Eigen::Map<Eigen::VectorXd> randE(rand.data(),rand.size());
    Eigen::VectorXd sampleE = meanHit + varHit * randE;
    std::vector<double> sample(sampleE.data(),sampleE.data() + sampleE.size());
    putPar(sample.begin());

    return true;
  }
  else{
    return false;
  }
}


void ModelBase::revert(){
  if(fitType == MLES){
    std::vector<double> sample(meanHit.data(),meanHit.data() + meanHit.size());
    putPar(sample.begin());
  }
  else{
    std::cout << "ModelBase::revert(): not implemented for fytType of "
	      << fitType << std::endl;
    throw(1);
  }
}
