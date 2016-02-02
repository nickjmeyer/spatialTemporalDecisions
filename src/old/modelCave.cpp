#include "modelCave.hpp"


void CaveModel::load(const SimData & sD,
		     const TrtData & tD,
		     const FixedData & fD,
		     const DynamicData & dD){
  mP.infProbsBase.zeros(sD.numInfected,sD.numNotInfec);
  mP.infProbsSep.zeros(sD.numInfected,sD.numNotInfec);

  arma::mat::iterator itBase;
  itBase = mP.infProbsBase.begin();
  
  int i,j,notNode;
  for(i=0; i<sD.numNotInfec; i++){
    notNode = sD.notInfec.at(i);
    for(j=0; j<sD.numInfected; j++,itBase++)
      (*itBase)=oneOnOne(notNode,sD.infected.at(j),sD,tD,fD,dD);
  }
  
  mP.setAll();
  
  mP.infProbs = arma::conv_to<std::vector<double> >
    ::from(1-arma::prod(mP.infProbsSep,0));
}



void CaveModel::infProbs(const SimData & sD,
			 const TrtData & tD,
			 const FixedData & fD,
			 const DynamicData & dD){
  mP.infProbs.clear();
  int i,j,node0;
  double prob;
  for(i=0; i<sD.numNotInfec; i++){
    node0 = sD.notInfec.at(i);
    prob=1.0;
    for(j=0; j<sD.numInfected; j++)
      prob *= 1/(1+std::exp(oneOnOne(node0,sD.infected.at(j),
				     sD,tD,fD,dD)));
    mP.infProbs.push_back(1-prob);
  }
}




void CaveModel::update(const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD){
  int i,j,k,node0,numNewInfec=sD.newInfec.size();
  double prob;
  std::vector<double> newInfProbs;
  for(i=0,j=0; i < sD.numNotInfec; ){ // purposely don't increment i or j here
    node0 = sD.notInfec.at(i);
    if(j == numNewInfec || node0 < sD.newInfec.at(j)){
      prob = 1 - mP.infProbs.at(i+j);
      for(k=0; k<numNewInfec; k++)
	prob*= 1/(1+std::exp(oneOnOne(node0,sD.newInfec.at(k),
				      sD,tD,fD,dD)));
      newInfProbs.push_back(1.0 - prob);
      i++;
    }
    else{
      j++;
    }
  }

  mP.infProbs = newInfProbs;
}



double CaveModel::oneOnOne(const int notNode,
			   const int infNode,
			   const SimData & sD,
			   const TrtData & tD,
			   const FixedData & fD,
			   const DynamicData & dD) const {
  double base = mP.intcp;

  base += mP.cave * fD.covar.at(notNode*fD.numCovar + 0);

  if(tD.p.at(notNode))
    base -= mP.trtPre;
  if(tD.a.at(infNode))
    base -= mP.trtAct;

  return base;
}



void CaveModel::fit(const SimData & sD, const TrtData & tD,
		    const FixedData & fD, const DynamicData & dD,
		    const int & useInit){
  if(useInit){
    fit(sD,tD,fD,dD,mP.getPar());
  }
  else{
    CaveParam mPInit;
    std::vector<double> par;
    par.push_back(-3.0);
    par.push_back(0.0);
    par.push_back(0.0);
    par.push_back(0.0);
  
    mPInit.putPar(par);
    fit(sD,tD,fD,dD,mPInit.getPar());
  }
}

void CaveModel::fit(const SimData & sD, const TrtData & tD,
		    const FixedData & fD, const DynamicData & dD,
		    const std::vector<double> & mPV){
  if(fitType == MLE){
    size_t iter=0;
    int status;

    gsl_vector *x,*ss;
    CaveParam mPInit;
    mPInit.putPar(mPV);
    std::vector<double> par = mPV;
    int i,dim=par.size();
    std::vector< std::vector<int> > history;
    history=sD.history;
    history.push_back(sD.status);
    CaveModelFitData dat(*this,mPInit,sD,fD,history);

    x = gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(x,i,par.at(i));
    ss=gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(ss,i,1.0);

    gsl_multimin_function minex_func;
    minex_func.n=dim;
    minex_func.f=&CaveModelFitObjFn;
    minex_func.params=&dat;

    const gsl_multimin_fminimizer_type *T=
      gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    s=gsl_multimin_fminimizer_alloc(T,dim);
    gsl_multimin_fminimizer_set(s,&minex_func,x,ss);

    double curSize;
    double size=.001;
  
    do{
      iter++;
      status=gsl_multimin_fminimizer_iterate(s);
      if(status)
	break;
      curSize=gsl_multimin_fminimizer_size(s);
      status=gsl_multimin_test_size(curSize,size);
    }while(status==GSL_CONTINUE && iter < 1000);

    for(i=0; i<dim; i++)
      par.at(i) = gsl_vector_get(s->x,i);
    mP.putPar(par);

    // if(sD.time <= fD.trtStart)
    //   mP.trtPre = mP.trtAct = 4.0;

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);

    load(sD,tD,fD,dD);
  }
  else if(fitType == MCMC){
    mcmc.load(sD.history,sD.status,fD);
    mcmc.sample(5000,1000,mPV);

    mcmc.samples.setMean();
    mP.putPar(mcmc.samples.getPar());

    load(sD,tD,fD,dD);
  }
  else{
    std::cout << "Not a valid Estimation" << std::endl;
    throw(1);
  }
}


CaveModelFitData
::CaveModelFitData(const CaveModel & m,
		   const CaveParam & mP,
		   const SimData & sD,
		   const FixedData & fD,
		   const std::vector<std::vector<int> > & history){
  this->m = m;
  this->mP = mP;
  this->fD = fD;
  this->history = history;
}

double CaveModelFitObjFn (const gsl_vector * x, void * params){
  CaveModelFitData * dat = static_cast<CaveModelFitData*> (params);
  double llike=0,prob,base;
  int i,j,t,time=dat->history.size(),dim=dat->mP.getPar().size();
  std::vector<double> par;

  for(i=0; i<dim; i++)
    par.push_back(gsl_vector_get(x,i));
  dat->mP.putPar(par);
  
  for(t=1; t<time; t++){
    for(i=0; i<dat->fD.numNodes; i++){
      if(dat->history.at(t-1).at(i) < 2){
	prob=1.0;
	for(j=0; j<dat->fD.numNodes; j++){
	  if(dat->history.at(t-1).at(j) >= 2){
	    base=dat->mP.intcp;
	    
	    base+=dat->mP.cave*dat->fD.caves.at(i);

	    if(dat->history.at(t-1).at(i) == 1)
	      base-=dat->mP.trtPre;
	    if(dat->history.at(t-1).at(j) == 3)
	      base-=dat->mP.trtAct;
	    
	    prob*=1.0/(1.0+std::exp(base));
	  }
	}
	if(dat->history.at(t).at(i)<2)
	  if(prob == 0)
	    llike+=-30;
	  else
	    llike+=std::log(prob);
	else
	  if(prob == 1)
	    llike+=-30;
	  else
	    llike+=std::log(1-prob);
      }
    }
  }
  
  if(!std::isfinite(llike))
    llike = std::numeric_limits<double>::lowest();
  
  return -llike;
}


