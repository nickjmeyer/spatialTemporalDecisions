#include "qOptimGSL.hpp"

QOptimTunePar::QOptimTunePar(){
  polReps = 2;
  radius = 50;

  jitter = .1;
  tol = .005;
  rate = 5;
  rateDecay = .975;
}

std::vector<double> QOptimTunePar::getPar() const{
  return std::vector<double> (0);
}

void QOptimTunePar::putPar(const std::vector<double> & par){
}



template class QOptim<System,RankAgent,GravityModel,GravityParam>;

template class QOptim<System,RankToyAgent,GravityModel,GravityParam>;


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
const int QOptim<System,Agent,Model,ModelParam>::numFeatures = 4;


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptim<System,Agent,Model,ModelParam>::
optim(System<Model,ModelParam> system,
      Agent<Model,ModelParam> & agent){
  
  std::vector<double> par=agent.tp.getPar();
  int i,sameRep=0,converged=0,numPar = par.size();
  std::vector<double> parJit(numPar);
  std::vector<double> parAdd(numPar);
  std::vector<double> parNew(numPar);

  
  double curr,prev = Qfn(system,agent);
  int iter=0;
  while(!converged){

    for(i=0; i<numPar; i++){
      parAdd.at(i) = tp.jitter*njm::rnorm01();
      parJit.at(i) = par.at(i) + parAdd.at(i);
    }

    agent.tp.putPar(par);
    prev = Qfn(system,agent);
    
    agent.tp.putPar(parJit);
    curr = Qfn(system,agent);

    njm::message("iter: " + njm::toString(iter++,"",4,0) +
    		 " || " + njm::toString(prev,"",6,4) + " - " +
    		 njm::toString(curr,"",6,4) + " -> " +
    		 njm::toString(par,", ",""));
    
    
    for(i=0; i<numPar; i++)
      parNew.at(i) = par.at(i) + tp.rate * (curr - prev) * parJit.at(i);

    njm::l2norm(parNew); // normalize
    
    if(njm::l2norm(parNew,par) < tp.tol){ // check if difference is small
      sameRep++;
      if(sameRep < 0)
	sameRep = 0;
      else if(sameRep == 3) // if small difference 3 times, it converged
	converged = 1;
    }

    tp.rate *= tp.rateDecay;
    par = parNew;
  }

  agent.tp.putPar(par); // assign optimized par to the agent

}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double QOptim<System,Agent,Model,ModelParam>::
bmRes(){
  return arma::as_scalar(R.t()*R + 2*R.t()*D*beta + beta.t()*D.t()*D*beta);
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
arma::colvec QOptim<System,Agent,Model,ModelParam>::
bmResG(){
  return 2*D.t()*R + 2*D.t()*D*beta;
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double QOptim<System,Agent,Model,ModelParam>::
spPen(){
  return arma::as_scalar(beta.t()*hmg*beta);
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
arma::colvec QOptim<System,Agent,Model,ModelParam>::
spPenG(){
  return 2*hmg*beta;
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptim<System,Agent,Model,ModelParam>::
setRD(const SimData & sD,
      const TrtData & tD,
      const FixedData & fD,
      const DynamicData & dD,
      const Model & m,
      ModelParam & mP,
      Agent<Model,ModelParam> & agent){
  SimData sDT;
  TrtData tDT;

  std::vector< std::vector<int> > h;
  h = sD.history;
  h.push_back(sD.status);

  // setup initial SimData
  sDT.time=0;
  sDT.status = h.at(0);
  sDT.timeInf.resize(fD.numNodes);
  std::fill(sDT.timeInf.begin(),sDT.timeInf.end(),0);
  int i;
  for(i=0; i<fD.numNodes; i++){
    if(sDT.status.at(i) < 2)
      sDT.notInfec.push_back(i);
    else{
      sDT.infected.push_back(i);
      sDT.newInfec.push_back(i);
      sDT.timeInf.at(i)++;
    }
  }
  sDT.numInfected = sDT.infected.size();
  sDT.numNotInfec = sDT.notInfec.size();

  tDT.a.resize(fD.numNodes);
  tDT.p.resize(fD.numNodes);
  std::fill(tDT.a.begin(),tDT.a.end(),0);
  std::fill(tDT.p.begin(),tDT.p.end(),0);

  arma::colvec y,dQ,ft,ft1;

  R.zeros(numFeatures*fD.numNodes);
  D.zeros(numFeatures*fD.numNodes,numFeatures*fD.numNodes);
  
  int t;
  for(t=0; t<(sD.time-1); t++){
    std::fill(tDT.a.begin(),tDT.a.end(),0);
    std::fill(tDT.p.begin(),tDT.p.end(),0);
    
    ft = getFeatures(sDT,tDT,fD,dD,m,mP);
    
    // next time point
    y.zeros(fD.numNodes);
    sDT.time++;
    sDT.history.push_back(sDT.status);
    sDT.status = h.at(t+1);
    sDT.notInfec.clear();
    sDT.infected.clear();
    sDT.newInfec.clear();
    for(i=0; i<fD.numNodes; i++){
      if(sDT.status.at(i) < 2)
	sDT.notInfec.push_back(i);
      else{
	sDT.infected.push_back(i);
	sDT.timeInf.at(i)++;
	if(h.at(t).at(i) < 2){
	  sDT.newInfec.push_back(i);
	  y(i) = 1;
	}
      }
    }
    sDT.numInfected = sDT.infected.size();
    sDT.numNotInfec = sDT.notInfec.size();

    if(sDT.time >= fD.trtStart){
      ft1.zeros(numFeatures*fD.numNodes);
      for(i=0; i<tp.polReps; i++){
    	agent.applyTrt(sDT,tDT,fD,dD,m,mP);
    	ft1 += getFeatures(sDT,tDT,fD,dD,m,mP);
      }
      ft1/=(double)tp.polReps;
    }
    else
      ft1 = getFeatures(sDT,tDT,fD,dD,m,mP);
    
    dQ = ft; // since our function is linear

    R += dQ * arma::sum(y);

    // discount factor is 1.0.  If it is changed < 1.0 must add an intercept
    D += dQ * (ft1 - ft).t();
  }
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptim<System,Agent,Model,ModelParam>::
setHMG(const SimData & sD, const FixedData & fD){
  int i,j,k,numNeigh,node0;
  arma::sp_mat H,G;

  hmg.zeros(fD.numNodes*numFeatures,fD.numNodes*numFeatures);

  std::vector< std::vector<int> > neighbors(fD.numNodes);
  numNeigh=0;
  for(i=0; i<fD.numNodes; i++)
    for(j=(i+1); j<fD.numNodes; j++)
      if(fD.dist.at(i*fD.numNodes + j) < tp.radius)
	neighbors.at(i).push_back(j);
    
  
  for(i=0; i<fD.numNodes; i++){
    numNeigh=neighbors.at(i).size();
    G.zeros(numFeatures*numNeigh,fD.numNodes*numFeatures);
    H.zeros(numFeatures*numNeigh,fD.numNodes*numFeatures);
    for(j=0; j<numNeigh; j++){
      node0=neighbors.at(i).at(j);
      for(k=0; k<numFeatures; k++){
	H(j*numFeatures+k,node0*numFeatures+k)=1;
	G(j*numFeatures+k,i*numFeatures+k)=1;
      }
    }

    if(numNeigh)
      hmg+= (H-G).t()*(H-G);
  }  
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
arma::colvec QOptim<System,Agent,Model,ModelParam>::
getFeatures(const SimData & sD,
	    const TrtData & tD,
	    const FixedData & fD,
	    const DynamicData & dD,
	    const Model & m,
	    ModelParam & mP){
  arma::mat infFeat,notFeat;
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);

  m.load(sD,tD,fD,dD,mP);

    int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;


  // feature 0
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,0).t();
  
  featNum++;
  
  // feature 1
  SystemLight<Model,ModelParam> s(sD,tD,fD,dD,m,mP);
  std::vector<int> newInfec;
  int k,numNewInfec,reps=2;
  for(i=0; i<reps; i++){
    s.reset();
    s.nextPoint();
    
    newInfec=s.sD.newInfec;
    s.nextPoint(1);
    
    newInfec.insert(newInfec.end(),s.sD.newInfec.begin(),s.sD.newInfec.end());
    
    numNewInfec = newInfec.size();
    for(j=0,itD0=newInfec.begin(); j<numNewInfec; j++,itD0++)
      for(k=0,itD1=sD.notInfec.begin(); k<sD.numNotInfec; k++,itD1++)
	if(*itD0 == *itD1)
	  notFeat(k,featNum)+=1.0/(double)reps;
  }

  infFeat.col(featNum) = (mP.infProbsSep * notFeat.col(featNum))
    /arma::sum(mP.infProbsSep,1);
  
  featNum++;

  
  // feature 2
  std::vector<double> infNoTrtLat,infNoTrtLong,notNoTrtLat,notNoTrtLong;
  int numInfNoTrt=0,numNotNoTrt=0;
  for(i=0; i<fD.numNodes; i++){
    if(sD.status.at(i) == 0){
      notNoTrtLat.push_back(fD.centroidsLat.at(i));
      notNoTrtLong.push_back(fD.centroidsLong.at(i));
      numNotNoTrt++;
    }
    else if(sD.status.at(i) == 2){
      infNoTrtLat.push_back(fD.centroidsLat.at(i));
      infNoTrtLong.push_back(fD.centroidsLong.at(i));
      numInfNoTrt++;
    }
  }

  
  for(i=0,itD0=sD.notInfec.begin(); i<sD.numNotInfec; i++,itD0++)
    notFeat(i,featNum) = halfPlaneDepth(fD.centroidsLong.at(*itD0),
					fD.centroidsLat.at(*itD0),
					numInfNoTrt,
					infNoTrtLong,
					infNoTrtLat);
  for(i=0,itD0=sD.infected.begin(); i<sD.numInfected; i++,itD0++)
    infFeat(i,featNum) = halfPlaneDepth(fD.centroidsLong.at(*itD0),
					fD.centroidsLat.at(*itD0),
					numNotNoTrt,
					notNoTrtLong,
					notNoTrtLat);

  featNum++;
  

  // feature 3
  std::vector<int>::const_iterator itD2,itD3;
  double totalDist,sigma=906.3819; // standard deviation of wnsD.txt
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      totalDist += std::exp(fD.dist.at((*itD0)*fD.numNodes + *itD1)/sigma);
    totalDist /= fD.numNodes*sigma;
    notFeat(i,featNum) = totalDist;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      totalDist += std::exp(fD.dist.at((*itD1)*fD.numNodes + *itD0)/sigma);
    totalDist /= fD.numNodes*sigma;
    infFeat(i,featNum) = totalDist;
  }

  featNum++;


  infFeat = infFeat.t();
  notFeat = notFeat.t();
  arma::mat::iterator notBeg,infBeg,notEnd,infEnd;
  notBeg = notFeat.begin();
  notEnd = notBeg + numFeatures;
  infBeg = infFeat.begin();
  infEnd = infBeg + numFeatures;

  std::vector<double> features;
  for(i=0,j=0; i<fD.numNodes; i++){
    if(j<sD.numInfected && sD.infected.at(j) == i){
      features.insert(features.end(),infBeg,infEnd);
      infBeg = infEnd;
      infEnd += numFeatures;
      j++;
    }
    else{
      features.insert(features.end(),notBeg,notEnd);
      notBeg = notEnd;
      notEnd += numFeatures;
    }
  }

  return arma::colvec(features);
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double QOptim<System,Agent,Model,ModelParam>::
Qfn(System<Model,ModelParam> s,
    Agent<Model,ModelParam> a){

  if((int)beta.n_elem != numFeatures*s.fD.numNodes)
    beta.zeros(numFeatures*s.fD.numNodes);
  
  QOptimData<System,Agent,Model,ModelParam> qD;
  qD.q.beta = beta;
  qD.q.tp = tp;
  qD.s = s;
  qD.a = a;
  
  size_t iter = 0;
  int status,dim=qD.q.beta.n_elem,i;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *w;
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  x=gsl_vector_alloc(dim);
  for(i=0; i<dim; i++)
    gsl_vector_set(x,i,qD.q.beta(i));
  
  my_func.n = dim;
  my_func.f = QOptimObjFn<System,Agent,Model,ModelParam>;
  my_func.df = QOptimObjG<System,Agent,Model,ModelParam>;
  my_func.fdf = QOptimObjFnG<System,Agent,Model,ModelParam>;
  my_func.params = &qD;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  w = gsl_multimin_fdfminimizer_alloc(T,dim);
  gsl_multimin_fdfminimizer_set(w,&my_func,x,0.1,0.1);

  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(w);

    if(status)
      break;

    if(iter < 10)
      status = GSL_CONTINUE;
    else
      status = gsl_multimin_test_gradient(w->gradient,1.0);


    // std::cout << "iter: " << std::setw(8) << iter << " -> " <<
    //   std::setw(8) << std::setprecision(4) << w->f << "\r" << std::flush;
  }
  while(status == GSL_CONTINUE && iter < 100);

  for(i=0; i<dim; i++)
    beta(i) = gsl_vector_get(w->x,i);

  gsl_multimin_fdfminimizer_free(w);
  gsl_vector_free(x);

  a.applyTrt(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);
  return(arma::as_scalar(getFeatures(s.sD,s.tD,s.fD,s.dD,
				     s.model,s.estParam).t()*beta));
}




template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double QOptimObjFn(const gsl_vector * x, void * param){
  QOptimData<System,Agent,Model,ModelParam>
    * qod = static_cast< QOptimData<System,Agent,Model,ModelParam> * > (param);
  int i,dim=qod->q.beta.n_elem;
  for(i=0; i<dim; i++)
    qod->q.beta(i) = gsl_vector_get(x,i);
  fixRandomSeed(1);
  qod->q.setRD(qod->s.sD,qod->s.tD,qod->s.fD,qod->s.dD,
	       qod->s.model,qod->s.estParam,qod->a);
  qod->q.setHMG(qod->s.sD,qod->s.fD);
  fixRandomSeed(0);
  return(qod->q.bmRes() + qod->q.spPen());
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptimObjG(const gsl_vector * x, void * param, gsl_vector * g){
  QOptimData<System,Agent,Model,ModelParam>
    * qod = static_cast< QOptimData<System,Agent,Model,ModelParam> * > (param);
  int i,dim=qod->q.beta.n_elem;
  for(i=0; i<dim; i++)
    qod->q.beta(i) = gsl_vector_get(x,i);
  fixRandomSeed(1);
  qod->q.setRD(qod->s.sD,qod->s.tD,qod->s.fD,qod->s.dD,
	       qod->s.model,qod->s.estParam,qod->a);
  qod->q.setHMG(qod->s.sD,qod->s.fD);
  fixRandomSeed(0);
  arma::colvec g_ = qod->q.bmResG() + qod->q.spPenG();
  for(i=0; i<dim; i++)
    gsl_vector_set(g,i,g_(i));
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void QOptimObjFnG(const gsl_vector * x, void * param,
		  double * f, gsl_vector * g){
  QOptimData<System,Agent,Model,ModelParam>
    * qod = static_cast< QOptimData<System,Agent,Model,ModelParam> * > (param);
  int i,dim=qod->q.beta.n_elem;
  for(i=0; i<dim; i++)
    qod->q.beta(i) = gsl_vector_get(x,i);
  fixRandomSeed(1);
  qod->q.setRD(qod->s.sD,qod->s.tD,qod->s.fD,qod->s.dD,
	       qod->s.model,qod->s.estParam,qod->a);
  qod->q.setHMG(qod->s.sD,qod->s.fD);
  fixRandomSeed(0);
  *f = qod->q.bmRes() + qod->q.spPen();
  arma::colvec g_ = qod->q.bmResG() + qod->q.spPenG();
  for(i=0; i<dim; i++)
    gsl_vector_set(g,i,g_(i));
}



