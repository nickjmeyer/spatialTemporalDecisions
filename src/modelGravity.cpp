#include "modelGravity.hpp"

static std::vector<ParamBase *> genPars(){
  std::vector<ParamBase *> pars;
  pars.push_back(new ParamIntercept);
  pars.push_back(new ParamBeta);
  pars.push_back(new ParamGravity);
  pars.push_back(new ParamTrt);
  return pars;
}

ModelGravity::ModelGravity(const FixedData & fD)
  : ModelBase(genPars(),fD){
}


ModelGravity::ModelGravity(const ModelGravity & m){
  int i, parsSize = m.pars.size();
  pars.clear();
  for(i = 0; i < parsSize; ++i)
    pars.push_back(m.pars.at(i)->clone());

  numPars = m.numPars;
  set = m.set;
  probs = m.probs;
  expitInfProbs = m.expitInfProbs;
  expitRevProbs = m.expitRevProbs;
  quick = m.quick;
  ready = m.ready;
  numInfected = m.numInfected;
  numNotInfec = m.numNotInfec;
  fitType = m.fitType;
  mcmc = m.mcmc;
}


ModelGravity & ModelGravity::operator=(const ModelGravity & m){
  if(this != & m){
    this->ModelGravity::~ModelGravity();
    new (this) ModelGravity(m);
  }
  return *this;
}



void ModelGravity::read(){
  std::vector<double> pars,add;
  njm::fromFile(add,njm::sett.srcExt("./GravityParam/intcp.txt"));
  pars.insert(pars.end(),add.begin(),add.end());
  
  njm::fromFile(add,njm::sett.srcExt("./GravityParam/beta.txt"));
  pars.insert(pars.end(),add.begin(),add.end());
  
  njm::fromFile(add,njm::sett.srcExt("./GravityParam/alpha.txt"));
  pars.insert(pars.end(),add.begin(),add.end());
  
  njm::fromFile(add,njm::sett.srcExt("./GravityParam/power.txt"));
  pars.insert(pars.end(),add.begin(),add.end());
  
  njm::fromFile(add,njm::sett.srcExt("./GravityParam/trtAct.txt"));
  pars.insert(pars.end(),add.begin(),add.end());
  
  njm::fromFile(add,njm::sett.srcExt("./GravityParam/trtPre.txt"));
  pars.insert(pars.end(),add.begin(),add.end());

  putPar(pars.begin());
}


void ModelGravity::save() const {
  std::vector<double> par;
  par = pars.at(0)->getPar();
  njm::toFile(njm::toString(par.at(0),"\n"),
	      njm::sett.srcExt("./GravityParam/intcp.txt"),
	      std::ios_base::out);

  par = pars.at(1)->getPar();
  njm::toFile(njm::toString(par,"\n",""),
	      njm::sett.srcExt("./GravityParam/beta.txt"),
	      std::ios_base::out);

  par = pars.at(2)->getPar();
  njm::toFile(njm::toString(par.at(0),"\n"),
	      njm::sett.srcExt("./GravityParam/alpha.txt"),
	      std::ios_base::out);
  
  njm::toFile(njm::toString(par.at(1),"\n"),
	      njm::sett.srcExt("./GravityParam/power.txt"),
	      std::ios_base::out);

  par = pars.at(3)->getPar();
  njm::toFile(njm::toString(par.at(0),"\n"),
	      njm::sett.srcExt("./GravityParam/trtAct.txt"),
	      std::ios_base::out);
  
  njm::toFile(njm::toString(par.at(1),"\n"),
	      njm::sett.srcExt("./GravityParam/trtPre.txt"),
	      std::ios_base::out);
}


double ModelGravity::tuneTrt(const FixedData & fD){
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

  std::vector<double> vals;
  vals = pars.at(0)->getPar();
  double base = vals.at(0);
  vals = pars.at(2)->getPar();
  double alpha = vals.at(0);
  double power = vals.at(1);
  base -= alpha * minDist/std::pow(avgCaves*avgCaves,power);

  return -(std::log(0.005) - base)/2.0;
}



void ModelGravity::fit(const SimData & sD, const TrtData & tD,
		       const FixedData & fD, const DynamicData & dD,
		       const int & useInit){
  if(useInit){
    fit(sD,tD,fD,dD,getPar());
  }
  else{
    std::vector<double> all;
    int i;
    for(i=0; i<(5+fD.numCovar); i++)
      all.push_back(0);
    all.at(0) = -3.0;
    fit(sD,tD,fD,dD,all);
  }
}

void ModelGravity::fit(const SimData & sD, const TrtData & tD,
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
    ModelGravityFitData dat(*this,all,fD,history);

    x = gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(x,i,all.at(i));
    ss=gsl_vector_alloc(dim);
    gsl_vector_set_all(ss,.5);

    gsl_multimin_function minex_func;
    minex_func.n=dim;
    minex_func.f=&modelGravityFitObjFn;
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
    } while(status==GSL_CONTINUE && iter < 1000);

    for(i=0; i<dim; i++)
      all.at(i) = gsl_vector_get(s->x,i);
    
    // njm::message(njm::toString(mP.getPar()," ","\n") +
    // 		 njm::toString(par," ","\n"));
    
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
    mcmc.load(sD.history,sD.status,fD);
    mcmc.sample(5000,1000,all);

    mcmc.samples.setMean();

    all = mcmc.samples.getPar();
    putPar(all.begin());

    setFill(sD,tD,fD,dD);
  }
  else{
    std::cout << "Not a valid Estimation of "
	      << fitType
	      << std::endl;
    throw(1);
  }
}


ModelGravityFitData
::ModelGravityFitData(const ModelGravity & m,
		      const std::vector<double> & all,
		      const FixedData & fD,
		      const std::vector<std::vector<int> > & history){
  this->m = m;
  this->m.putPar(all.begin());
  this->fD = fD;
  this->history = history;
}

double modelGravityFitObjFn (const gsl_vector * x, void * params){
  ModelGravityFitData * dat = static_cast<ModelGravityFitData*> (params);
  double llike=0,prob,base,caveTerm;
  int i,j,k,t,time=dat->history.size(),dim=dat->m.getPar().size();
  std::vector<double> par;

  for(i=0; i<dim; i++)
    par.push_back(gsl_vector_get(x,i));
  dat->m.putPar(par.begin());

  std::vector<double>::const_iterator it = par.begin();
  double intcp = *it;
  ++it;
  std::vector<double> beta;
  for(i = 0; i < dat->fD.numCovar; ++i,++it)
    beta.push_back(*it);
  double alpha = *it;
  ++it;
  double power = *it;
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
	    for(k=0; k<dat->fD.numCovar; k++)
	      base+=beta.at(k)*dat->fD.covar.at(i*dat->fD.numCovar+k);
	    caveTerm=dat->fD.dist.at(i*dat->fD.numNodes+j);
	    caveTerm/=std::pow(dat->fD.caves.at(i)*dat->fD.caves.at(j),
			       power);
	    base-=alpha*caveTerm;
	    if(dat->history.at(t-1).at(i) == 1)
	      base-=trtPre;
	    if(dat->history.at(t-1).at(j) == 3)
	      base-=trtAct;
	    prob*=1/(1+std::exp(base));
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



