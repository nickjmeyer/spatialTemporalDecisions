#include "modelRad.hpp"


static std::vector<ParamBase *> genPars(){
  std::vector<ParamBase *> pars;
  pars.push_back(new ParamIntercept);
  pars.push_back(new ParamBeta);
  pars.push_back(new ParamGDist);
  pars.push_back(new ParamRad);
  pars.push_back(new ParamTrt);
  return pars;
}

ModelRad::ModelRad(const FixedData & fD)
  : ModelBase("Rad",genPars(),fD){
}


ModelRad::ModelRad(const ModelRad & m){
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
}


ModelRad & ModelRad::operator=(const ModelRad & m){
  if(this != & m){
    this->ModelRad::~ModelRad();
    new (this) ModelRad(m);
  }
  return *this;
}



void ModelRad::fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit){
  if(useInit){
    fit(sD,tD,fD,dD,getPar());
  }
  else{
    std::vector<double> all;
    all.push_back(-3.0); // intcp
    for(int i = 0; i < fD.numCovar; ++i)
      all.push_back(0.0);
    all.push_back(0.0); // GDist
    all.push_back(0.0); // Rad
    all.push_back(0.0); // trt
    all.push_back(0.0); // trt
  
    fit(sD,tD,fD,dD,all);
  }
}

void ModelRad::fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all){
  if(fitType == MLE || fitType == MLES){
    size_t iter=0;
    int status;

    gsl_vector *x,*ss;
    int i,dim=all.size();
    std::vector< std::vector<int> > history;
    history=sD.history;
    history.push_back(sD.status);
    ModelRadFitData dat(*this,all,fD,history);

    x = gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(x,i,all.at(i));
    ss=gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(ss,i,1.0);

    gsl_multimin_function minex_func;
    minex_func.n=dim;
    minex_func.f=&modelRadFitObjFn;
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
      all.at(i) = gsl_vector_get(s->x,i);

    if(sD.time <= fD.trtStart){
      std::vector<double>::iterator it;
      it = all.end();
      --it;
      *it = fD.priorTrtMean;
      --it;
      *it = fD.priorTrtMean;
    }

    putPar(all.begin());
    

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);

    if(fitType == MLES)
      setFisher(sD,tD,fD,dD);

    setFill(sD,tD,fD,dD);
  }
  else if(fitType == MCMC){
    std::cout << "Error: ModelRad::fit(): MCMC not setup"
	      << std::endl;
    throw(1);
  }
  else{
    std::cout << "Not a valid Estimation" << std::endl;
    throw(1);
  }
}


ModelRadFitData
::ModelRadFitData(const ModelRad & m,
		  const std::vector<double> & all,
		  const FixedData & fD,
		  const std::vector<std::vector<int> > & history){
  this->m = m;
  this->m.putPar(all.begin());
  this->fD = fD;
  this->history = history;
}

double modelRadFitObjFn (const gsl_vector * x, void * params){
  ModelRadFitData * dat = static_cast<ModelRadFitData*> (params);
  double llike=0,prob,base;
  int i,j,t,time=dat->history.size(),dim=dat->m.getPar().size();
  std::vector<double> par;

  for(i=0; i<dim; i++)
    par.push_back(gsl_vector_get(x,i));
  dat->m.putPar(par.begin());

  std::vector<double>::const_iterator it = par.begin();
  double intcp = *it;
  ++it;
  double alpha = *it;
  ++it;
  double rad = *it;
  ++it;
  double trtAct = *it;
  ++it;
  double trtPre = *it;

  
  for(t=1; t<time; t++){
    for(i=0; i<dat->fD.numNodes; i++){
      if(dat->history.at(t-1).at(i) < 2){
	prob=1.0;
	for(j=0; j<dat->fD.numNodes; j++){
	  if(dat->history.at(t-1).at(j) >= 2){
	    base=intcp;
	    
	    base-=alpha*dat->fD.gDist.at(i*dat->fD.numNodes+j);

	    double n_i,n_j,cm_ij;
	    n_i = dat->fD.caves[i];
	    n_j = dat->fD.caves[j];
	    cm_ij = dat->fD.cm[i*dat->fD.numNodes + j];

	    base += rad * n_i*n_i*n_j/(n_i*cm_ij*(n_i+n_j+cm_ij));

	    if(dat->history.at(t-1).at(i) == 1)
	      base-=trtPre;
	    if(dat->history.at(t-1).at(j) == 3)
	      base-=trtAct;
	    
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
  return -llike;
}


