#include "modelGravityTimeInfExpCaves.hpp"


double
GravityTimeInfExpCavesModel::tuneTrt(const FixedData & fD,
				     const GravityTimeInfExpCavesParam & gP){
  int i,j;
  double avgCaves = 0.0;
  for(i = 0; i < fD.numNodes; i++)
    avgCaves += fD.caves.at(i);
  avgCaves /= double(fD.numNodes);

  double minDist = std::numeric_limits<double>::max();
  for(i = 0; i < fD.numNodes; i++)
    for(j = (i+1); j < fD.numNodes; j++)
      if(minDist > fD.dist.at(i*fD.numNodes + j))
	minDist = fD.dist.at(i*fD.numNodes + j);

  double base = gP.intcp;
  base -= gP.alpha * minDist/std::pow(avgCaves*avgCaves,gP.power);

  return -(std::log(0.005) - base)/2.0;
}


double
GravityTimeInfExpCavesModel::oneOnOne(const int notNode,
				      const int infNode,
				      const SimData & sD,
				      const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD,
				      const
				      GravityTimeInfExpCavesParam & gP) const{
  double base = gP.intcp;
  int len = gP.beta.size();
  for(int i=0; i<len; i++)
    base += gP.beta.at(i)*fD.covar.at(notNode*len + i);

  base -= gP.alpha * fD.dist.at(notNode*fD.numNodes + infNode)/
    std::pow(fD.caves.at(notNode)*fD.caves.at(infNode),gP.power);

  double caveTime = (sD.timeInf.at(infNode)-1.0)/fD.propCaves.at(infNode);
  base += gP.xi * (std::exp(caveTime)-1.0);
  
  if(tD.p.at(notNode))
    base -= gP.trtPre;
  if(tD.a.at(infNode))
    base -= gP.trtAct;

  return base;
}



void
GravityTimeInfExpCavesModel::fit(const SimData & sD, const TrtData & tD,
				 const FixedData & fD, const DynamicData & dD,
				 GravityTimeInfExpCavesParam & mP){
  GravityTimeInfExpCavesParam mPInit;
  std::vector<double> par;
  int i;
  for(i=0; i<(6+fD.numCovar); i++)
    par.push_back(0);
  mPInit.putPar(par);
  mPInit.intcp=-3.0;
  fit(sD,tD,fD,dD,mP,mPInit);
}

void GravityTimeInfExpCavesModel::fit(const SimData & sD, const TrtData & tD,
				      const FixedData & fD,
				      const DynamicData & dD,
				      GravityTimeInfExpCavesParam & mP,
				      const GravityTimeInfExpCavesParam mPInit){

  if(fitType == MLE){
  
    size_t iter=0;
    int status;

    gsl_vector *x,*ss;
    std::vector<double> par = mPInit.getPar();
    int i,dim=par.size();
    std::vector< std::vector<int> > history;
    history=sD.history;
    history.push_back(sD.status);
    GravityTimeInfExpCavesModelFitData dat(*this,mPInit,sD,fD,history);

    x = gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(x,i,par.at(i));
    ss=gsl_vector_alloc(dim);
    gsl_vector_set_all(ss,.5);

    gsl_multimin_function minex_func;
    minex_func.n=dim;
    minex_func.f=&gravityTimeInfExpCavesModelFitObjFn;
    minex_func.params=&dat;

    const gsl_multimin_fminimizer_type *T=
      gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    s=gsl_multimin_fminimizer_alloc(T,dim);
    gsl_multimin_fminimizer_set(s,&minex_func,x,ss);

    double curSize;
    double size=0.001;
  
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

    if(sD.time <= fD.trtStart)
      mP.trtPre = mP.trtAct = fD.priorTrtMean;

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);

    load(sD,tD,fD,dD,mP);
    
  }
  else if(fitType == MCMC){
    mcmc.load(sD.history,sD.status,fD);
    mcmc.sample(5000,1000,mPInit.getPar());

    mcmc.samples.setMean();
    mP.putPar(mcmc.samples.getPar());

    load(sD,tD,fD,dD,mP);
  }
  else{
    std::cout << "Not a valid Estimation" << std::endl;
    throw(1);
  }
}


GravityTimeInfExpCavesModelFitData
::GravityTimeInfExpCavesModelFitData(const GravityTimeInfExpCavesModel & m,
				     const GravityTimeInfExpCavesParam & mP,
				     const SimData & sD,
				     const FixedData & fD,
				     const
				     std::vector<std::vector<int> > & history){
  this->m = m;
  this->mP = mP;
  this->sD = sD;
  this->fD = fD;
  this->history = history;

  this->timeInf.clear();
  std::vector<int> timeInf(fD.numNodes,0);
  int i,j;
  for(i = 0; i < (int)history.size(); ++i){
    for(j = 0; j < fD.numNodes; ++j){
      if(history.at(i).at(j) >= 2)
	++timeInf.at(j);
    }
    this->timeInf.push_back(timeInf);
  }
}

double
gravityTimeInfExpCavesModelFitObjFn (const gsl_vector * x, void * params){
  GravityTimeInfExpCavesModelFitData * dat =
    static_cast<GravityTimeInfExpCavesModelFitData*> (params);
  double llike=0,prob,base,caveTerm;
  int i,j,k,t,time=dat->history.size(),dim=dat->mP.getPar().size();
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
	    for(k=0; k<dat->fD.numCovar; k++)
	      base+=dat->mP.beta.at(k)*dat->fD.covar.at(i*dat->fD.numCovar+k);
	    caveTerm=dat->fD.dist.at(i*dat->fD.numNodes+j);
	    caveTerm/=std::pow(dat->fD.caves.at(i)*dat->fD.caves.at(j),
			       dat->mP.power);
	    base-=dat->mP.alpha*caveTerm;

	    base+=dat->mP.xi*(std::exp((dat->timeInf.at(t-1).at(j) - 1.0)/
				       dat->fD.propCaves.at(j))-1.0);
	    
	    if(dat->history.at(t-1).at(i) == 1)
	      base-=dat->mP.trtPre;
	    if(dat->history.at(t-1).at(j) == 3)
	      base-=dat->mP.trtAct;
	    
	    prob*=1/(1+std::exp(base));
	  }
	}
	if(dat->history.at(t).at(i)<2)
	  if(prob == 0)
	    llike+=-30;
	  else
	    llike+=std::exp(prob);
	else
	  if(prob == 1)
	    llike+=-30;
	  else
	    llike+=std::exp(1-prob);
      }
    }
  }
  return -llike;
}



