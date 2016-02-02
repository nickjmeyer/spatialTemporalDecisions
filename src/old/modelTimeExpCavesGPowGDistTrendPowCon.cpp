#include "modelTimeExpCavesGPowGDistTrendPowCon.hpp"

static std::vector<ParamBase *> genPars(){
  std::vector<ParamBase *> pars;
  pars.push_back(new ParamIntercept);
  pars.push_back(new ParamBeta);
  pars.push_back(new ParamGravPowGDist);
  pars.push_back(new ParamTrendPowCon);
  pars.push_back(new ParamTimeExpCaves);
  pars.push_back(new ParamTrt);
  return pars;
}


ModelTimeExpCavesGPowGDistTrendPowCon
::ModelTimeExpCavesGPowGDistTrendPowCon(const FixedData & fD)
  : ModelBase("TimeExpCavesGPowGDistTrendPowCon",genPars(),fD){
}


ModelTimeExpCavesGPowGDistTrendPowCon
::ModelTimeExpCavesGPowGDistTrendPowCon
(const ModelTimeExpCavesGPowGDistTrendPowCon & m){

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
  fixSample = m.fixSample;
}


ModelTimeExpCavesGPowGDistTrendPowCon &
ModelTimeExpCavesGPowGDistTrendPowCon
::operator=(const ModelTimeExpCavesGPowGDistTrendPowCon & m){
  if(this != & m){
    this->ModelTimeExpCavesGPowGDistTrendPowCon
      ::~ModelTimeExpCavesGPowGDistTrendPowCon();
    new (this) ModelTimeExpCavesGPowGDistTrendPowCon(m);
  }
  return *this;
}



void
ModelTimeExpCavesGPowGDistTrendPowCon::fit(const SimData & sD,
					   const TrtData & tD,
					   const FixedData & fD,
					   const DynamicData & dD,
					   const int & useInit){
  if(useInit){
    fit(sD,tD,fD,dD,getPar());
  }
  else{
    std::vector<double> all;
    unsigned int i;
    for(i=0; i<numPars; i++)
      all.push_back(0);
    all.at(0) = -3.0;
    fit(sD,tD,fD,dD,all);
  }
}

void ModelTimeExpCavesGPowGDistTrendPowCon::fit(const SimData & sD,
					    const TrtData & tD,
					    const FixedData & fD,
					    const DynamicData & dD,
					    std::vector<double> all){

  if(fitType == MLE || fitType == MLES){

    size_t iter=0;
    int status;

    gsl_vector *x,*ss;
    int i,dim=all.size();
    std::vector< std::vector<int> > history;
    history=sD.history;
    history.push_back(sD.status);
    ModelTimeExpCavesGPowGDistTrendPowConFitData dat(*this,all,fD,history);

    x = gsl_vector_alloc(dim);
    for(i=0; i<dim; i++)
      gsl_vector_set(x,i,all.at(i));
    ss=gsl_vector_alloc(dim);
    gsl_vector_set_all(ss,.5);

    gsl_multimin_function minex_func;
    minex_func.n=dim;
    minex_func.f=&modelTimeExpCavesGPowGDistTrendPowConFitObjFn;
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
      all.at(i) = gsl_vector_get(s->x,i);

    if(sD.time <= fD.trtStart){
      std::vector<std::string> parNames = {"trtAct","trtPre"};
      setPar(parNames,fD.priorTrtMean);
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
    std::cout << "Error: ModelTimeExpCavesGPowGDistTrendPowCon"
	      << "::fit(): MCMC not setup"
	      << std::endl;
    throw(1);
  }
  else{
    std::cout << "Not a valid Estimation" << std::endl;
    throw(1);
  }
}


ModelTimeExpCavesGPowGDistTrendPowConFitData
::ModelTimeExpCavesGPowGDistTrendPowConFitData
(const ModelTimeExpCavesGPowGDistTrendPowCon & m,
 const std::vector<double> & all,
 const FixedData & fD,
 const
 std::vector<std::vector<int> > & history){

  this->m = m;
  this->m.putPar(all.begin());
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
modelTimeExpCavesGPowGDistTrendPowConFitObjFn(const gsl_vector * x,
					      void * params){
  ModelTimeExpCavesGPowGDistTrendPowConFitData * dat =
    static_cast<ModelTimeExpCavesGPowGDistTrendPowConFitData*> (params);
  double llike=0,prob,base,caveTerm;
  int i,j,k,t,time=dat->history.size(),dim=dat->m.getPar().size();

  std::vector<double> par;
  for(i=0; i<dim; i++)
    par.push_back(gsl_vector_get(x,i));

  std::vector<double>::const_iterator it = par.begin();

  double intcp = *it++;
  std::vector<double> beta;
  for(i = 0; i < dat->fD.numCovar; ++i)
    beta.push_back(*it++);
  double alpha = *it++;
  double power = *it++;
  double gPow = *it++;
  double trend = *it++;
  double trendPowCon = *it++;
  double xi = *it++;
  double trtAct = *it++;
  double trtPre = *it++;

  for(t=1; t<time; t++){
    for(i=0; i<dat->fD.numNodes; i++){
      if(dat->history.at(t-1).at(i) < 2){
	prob=1.0;
	for(j=0; j<dat->fD.numNodes; j++){
	  if(dat->history.at(t-1).at(j) >= 2){
	    base=intcp;
	    for(k=0; k<dat->fD.numCovar; k++)
	      base+=beta.at(k)*dat->fD.covar.at(i*dat->fD.numCovar+k);
	    caveTerm=std::pow(dat->fD.gDist.at(i*dat->fD.numNodes+j),gPow);
	    caveTerm/=std::pow(dat->fD.caves.at(i)*dat->fD.caves.at(j),
			       power);
	    base-=alpha*caveTerm;

	    base+=xi*(std::exp((dat->timeInf.at(t-1).at(j) - 1.0)/
			       dat->fD.propCaves.at(j))-1.0);

	    base+=trend*std::pow(double(t-1),-std::exp(trendPowCon));
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
	    llike+=std::exp(prob);
	else
	  if(prob == 1)
	    llike+=-30;
	  else
	    llike+=std::exp(1-prob);
      }
    }
  }

  if(!std::isfinite(llike))
    llike = std::numeric_limits<double>::lowest();

  return -llike;
}
