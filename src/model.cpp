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


void ModelBase::setFixSample(const int & fix){
  fixSample = fix;
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
  // njm::message("setFill");
  int i,numPars = pars.size();
  probs = std::vector<double>(fD.numNodes*fD.numNodes,0.0);
  for(i = 0; i < numPars; ++i){
    pars[i]->setFill(probs,sD,tD,fD,dD);
    // njm::message(njm::toString(i,"",0,0) + ": " +
    // 		 njm::toString(std::accumulate(probs.begin(),probs.end(),0.0),
    // 			       ""));
  }
  set = 1;
  ready = 0;
}


void ModelBase::modFill(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  if(set == 1){
    // njm::message("modFill");
    int i,numPars = pars.size();
    for(i = 0; i < numPars; ++i){
      pars[i]->modFill(probs,sD,tD,fD,dD);
      // njm::message(njm::toString(i,"",0,0) + ": " +
      // 		   njm::toString(std::accumulate(probs.begin(),
      // 						 probs.end(),0.0),
      // 				 ""));
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
  ready = 1;
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
  bool found = false;
  for(i = 0; i < numPars; ++i){
    found |= pars[i]->setPar(name,val);
  }

  if(!found){
    std::cout << "name " + name + " not found in setPar()" << std::endl;
    throw(1);
  }
}


void ModelBase::setPar(const std::vector<std::string> & name,
		       const double & val){
  int i,numPars = pars.size();
  bool found = false;
  for(i = 0; i < numPars; ++i)
    found |= pars[i]->setPar(name,val);

  if(!found){
    std::cout << "name not found in setPar()" << std::endl;
    throw(1);
  }
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

  // std::cout << "p: " << njm::toString(p) << std::endl;
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

  // std::cout << totLen << std::endl;
  // std::cout << p[9*totLen + 9] << ", "
  // 	    << p[9*totLen + 10] << ", " << p[9*totLen + 11] << std::endl;
  // std::cout << p[10*totLen + 9] << ", "
  // 	    << p[10*totLen + 10] << ", " << p[10*totLen + 11] << std::endl;
  // std::cout << p[11*totLen + 9] << ", "
  // 	    << p[11*totLen + 10] << ", " << p[11*totLen + 11] << std::endl;

  return p;
}


void ModelBase::setFisher(const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
  std::cout << "=========== fisher ============" << std::endl;
  fisher.resize(numPars * numPars);
  std::fill(fisher.begin(),fisher.end(),0);

  std::cout << njm::toString(getPar()) << std::endl;
  std::cout << "power: " << getPar({"power"})[0] << std::endl;

  std::vector<std::vector<int> > hist = sD.history;
  hist.push_back(sD.status);
  std::vector<DataBundle> db = historyToData(hist);

  int t,iN,nN;
  unsigned int pi,pj;
  for(t = 0; t < sD.time; ++t){
    SimData sDi= std::get<0>(db[t]);
    TrtData tDi= std::get<1>(db[t]);
    DynamicData dDi = std::get<2>(db[t]);

    std::cout << "t = " << t << std::endl;

    std::cout << "fill" << std::endl;
    setFill(sDi,tDi,fD,dDi);
    std::cout << "end fill" << std::endl;
    setQuick(sDi,tDi,fD,dDi);
    infProbs(sDi,tDi,fD,dDi);

    if(t == 0){
      std::cout << "infected: " << std::endl;
      std::cout << sDi.numInfected << std::endl;
      std::cout << njm::toString(sDi.infected) << std::endl;

      std::cout << "quick: " << std::endl;
      std::cout << njm::toString(quick) << std::endl;

      std::cout << "infProbs: " << std::endl;
      std::cout << njm::toString(expitInfProbs) << std::endl;

      std::cout << "neighbors of " << sDi.infected.at(0) << ": " << std::endl;
      for(int i = 0; i < fD.numNodes; ++i){
	if(fD.network.at(i*fD.numNodes + sDi.infected.at(0)) == 1){
	  std::cout << i << ": " << probs[i*fD.numNodes + sDi.infected.at(0)]
		    << " ++ " << fD.gDist[i*fD.numNodes + sDi.infected.at(0)]
		    << std::endl;
	}
      }
      std::cout << std::endl;

      // std::cout << "probs: " << std::endl;
      // std::cout << njm::toString(probs) << std::endl;
    }

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

	// double quickInd = std::max(1e-10,quick[ind]);
	double quickInd = 1.0 - quick[ind];
	for(pi = 0; pi < numPars; ++pi){
	  sqr[pi] += quickInd*p[pi];

	  int piInd = pi*numPars;
	  for(pj = pi; pj < numPars; ++pj){
	    dbl[piInd + pj] += quickInd*p2[piInd + pj];
	    dbl[piInd + pj] += quickInd*(1.0-quickInd)*p[pj]*p[pi];
	    if(pj != pi){
	      dbl[pj*numPars + pi] += quickInd*p2[pj*numPars + pi];
	      dbl[pj*numPars + pi] += quickInd*(1.0-quickInd)*p[pi]*p[pj];
	    }
	  }
	}
      }

      int next = (hist[t+1][sDi.notInfec[nN]] < 2 ? 0 : 1);

      for(pi = 0; pi < numPars; ++pi){
	for(pj = pi; pj < numPars; ++pj){
	  // double prob = std::max(1e-10,expitInfProbs[nN]);
	  double prob = std::max(1e-10,expitInfProbs[nN]);

	  if(prob > 0.0){
	    fisher[pi*numPars + pj] +=
	      (double(next)/prob - 1.0)*dbl[pi*numPars + pj];

	    if(pj != pi){
	      fisher[pj*numPars + pi] +=
		(double(next)/prob - 1.0)*dbl[pj*numPars + pi];
	    }
	  }

	  if((prob*prob) > 0.0){
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

  for(int i = 0; i < int(I.rows()); ++i){
    std::cout << "I row " << i << ": " << I.row(i) << std::endl;
  }

  if(ldlt.info() != Eigen::Success){
    static int badCnt = 0;
    std::cout << "Error: ModelBase::getFisher(): eigen solver failed"
	      << std::endl;
    njm::toFile(fisher,njm::sett.datExt("badFisher" +
					njm::toString(badCnt++,"",0,0)));
    throw(1);
  }

  else{
    std::cout << "here" << std::endl;
    std::cout << I(0,0) << std::endl;
    Eigen::Transpositions<Eigen::Dynamic,Eigen::Dynamic> P;
    P = ldlt.transpositionsP();
    Eigen::MatrixXd L = ldlt.matrixL();

    for(int i = 0; i < int(L.rows()); ++i){
      std::cout << "L row " << i << ": " << L.row(i) << std::endl;
    }

    // std::cout << "L: " << std::endl
    // 	      << L << std::endl
    // 	      << std::endl;
    L.transposeInPlace();
    Eigen::VectorXd D = ldlt.vectorD();

    std::cout << "D: "  << std::endl << D << std::endl;

    throw(1);

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


bool ModelBase::sample(const bool force){
  if(fitType == MLES && (!fixSample || force)){
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




void ModelBase::estimateMle(const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){
  std::vector<double> startingVals = this->getPar();
  estimateMle(startingVals,sD,tD,fD,dD);
}


void ModelBase::estimateMle(const std::vector<double> startingVals,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD){

  if(fitType != MLE && fitType != MLES){
    std::cout << "ModelBase::estimateMLE is called for non MLE or MLES type"
	      << std::endl;
    throw(1);
  }

  ModelBaseFitObj fitObj(this,sD,tD,fD,dD);

  const gsl_multimin_fdfminimizer_type * T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector * x;
  x = gsl_vector_alloc(this->numPars);
  int pi;
  for(pi = 0; pi < int(this->numPars); ++pi){
    gsl_vector_set(x,pi,startingVals.at(pi));
  }

  gsl_multimin_function_fdf my_func;
  my_func.n = this->numPars;
  my_func.f = objFn;
  my_func.df = objFnGrad;
  my_func.fdf = objFnBoth;
  my_func.params = &fitObj;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc(T,this->numPars);

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.01,1e-4);

  int iter = 0;
  int status;
  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if(status)
      break;

    status = gsl_multimin_test_gradient(s->gradient,1e-3);


    if(status == GSL_SUCCESS){
      std::cout << "Min found" << std::endl;
    }
  }while(status == GSL_CONTINUE && iter < 100);


  std::vector<double> mle;
  for(pi = 0; pi < int(numPars); ++pi){
    mle.push_back(gsl_vector_get(s->x,pi));
  }

  std::cout << "mle: " << njm::toString(mle) << std::endl;
  std::cout << "par: " << njm::toString(getPar()) << std::endl;
  this->putPar(mle.begin());


  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
}



double ModelBase::logll(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){

  std::vector<std::vector<int> > hist = sD.history;
  hist.push_back(sD.status);
  std::vector<DataBundle> db = historyToData(hist);

  double logllVal = 0.0;

  int t,nN;
  // loop over time points
  for(t = 0; t < sD.time; ++t){
    SimData sDi = std::get<0>(db[t]);
    TrtData tDi = std::get<1>(db[t]);
    DynamicData dDi = std::get<2>(db[t]);

    setFill(sDi,tDi,fD,dDi);
    infProbs(sDi,tDi,fD,dDi);

    if(int(expitInfProbs.size()) != sDi.numNotInfec){
      std::cout << "ModelBase::logll(): length of expitInfProbs is not same as"
		<< " number of uninfected nodes at time t"
		<< std::endl;
      throw(1);
    }

    // loop over uninfected nodes at time t
    for(nN = 0; nN < sDi.numNotInfec; ++nN){
      double prob = expitInfProbs.at(nN);
      int next = (hist[t+1][sDi.notInfec[nN]] < 2) ? 0 : 1;
      if(next == 1){
	if(prob < 1e-44)
	  logllVal += -100.0;
	else
	  logllVal += std::log(prob);
      }
      else{
	if((1.0-prob) < 1e-44)
	  logllVal += -100.0;
	else
	  logllVal += std::log(1.0 - prob);
      }
    }
  }
  return logllVal;
}


std::vector<double> ModelBase::logllGrad(const SimData & sD,
					 const TrtData & tD,
					 const FixedData & fD,
					 const DynamicData & dD){
  std::vector<std::vector<int> > hist = sD.history;
  hist.push_back(sD.status);
  std::vector<DataBundle> db = historyToData(hist);

  std::vector<double> logllGradVal(numPars,0.0);
  std::fill(logllGradVal.begin(),logllGradVal.end(),0.0);

  int t,nN,iN,pi;
  // loop over time points
  for(t = 0; t < sD.time; ++t){
    SimData sDi = std::get<0>(db[t]);
    TrtData tDi = std::get<1>(db[t]);
    DynamicData dDi = std::get<2>(db[t]);

    setFill(sDi,tDi,fD,dDi);
    setQuick(sDi,tDi,fD,dDi);
    infProbs(sDi,tDi,fD,dDi);

    if(int(expitInfProbs.size()) != sDi.numNotInfec){
      std::cout << "ModelBase::logll(): length of expitInfProbs is not same as"
		<< " number of uninfected nodes at time t"
		<< std::endl;
      throw(1);
    }

    // loop over uninfected nodes at time t
    for(nN = 0; nN < sDi.numNotInfec; ++nN){
      double prob = expitInfProbs.at(nN);
      int next = (hist[t+1][sDi.notInfec[nN]] < 2) ? 0 : 1;

      if(prob > 0.0){
	double beg = double(next)/prob - 1.0;
	for(iN = 0; iN < sDi.numInfected; ++iN){
	  // quick stores probabilty of not infefcting need to take
	  // (1.0 - quick) to get probablity of infecting

	  double nNInfByiN = 1.0 - quick[nN*sDi.numInfected + iN];
	  std::vector<double> grad = partial(sDi.notInfec[nN],
					     sDi.infected[iN],
					     sDi,tDi,fD,dDi);
	  for(pi=0; pi < int(numPars); ++pi){
	    logllGradVal.at(pi) += beg * nNInfByiN * grad.at(pi);
	  }
	}
      }
    }
  }
  return logllGradVal;
}



ModelBaseFitObj::ModelBaseFitObj(ModelBase * const mb,
				 const SimData sD,
				 const TrtData tD,
				 const FixedData fD,
				 const DynamicData dD){
  this->mb = mb;
  this->sD = sD;
  this->tD = tD;
  this->fD = fD;
  this->dD = dD;
}


double objFn(const gsl_vector * x, void * params){
  ModelBaseFitObj * fitObj = static_cast<ModelBaseFitObj*>(params);
  std::vector<double> par;
  int pi;
  for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
    par.push_back(gsl_vector_get(x,pi));
  }

  fitObj->mb->putPar(par.begin());

  // return negative since GSL minimizes the function
  double ll = fitObj->mb->logll(fitObj->sD,fitObj->tD,fitObj->fD,fitObj->dD);
  return - ll;
}

void objFnGrad(const gsl_vector * x, void * params, gsl_vector * g){
  ModelBaseFitObj * fitObj = static_cast<ModelBaseFitObj*>(params);
  std::vector<double> par;
  int pi;
  for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
    par.push_back(gsl_vector_get(x,pi));
  }

  fitObj->mb->putPar(par.begin());

  std::vector<double> llGrad = fitObj->mb->logllGrad(fitObj->sD,fitObj->tD,
						     fitObj->fD,fitObj->dD);
  for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
    // assign the negative of the gradient value
    // GSL minimizes the function, need to adjust the gradient too
    gsl_vector_set(g,pi,-llGrad.at(pi));
  }

}

void objFnBoth(const gsl_vector * x, void * params, double * f, gsl_vector * g){
  *f = objFn(x,params);
  objFnGrad(x,params,g);
}
