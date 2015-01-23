#include "m2NmOptimOld.hpp"

M2OldNmEvalTunePar::M2OldNmEvalTunePar(){
  polReps = 20;
  // radius = 50;
  valReps = 100;
  numNeigh = 5;
  
  gamma = .95;
  lambda = 450000.0;
  
  jitter = .1;
  tol = .005;
  rate = 5;
  rateDecay = .975;
}

std::vector<double> M2OldNmEvalTunePar::getPar() const{
  return std::vector<double> (0);
}

void M2OldNmEvalTunePar::putPar(const std::vector<double> & par){
}



template class M2OldNmOptim<System,RankAgent,GravityModel,GravityParam>;

template class M2OldNmOptim<System,RankToyAgent,GravityModel,GravityParam>;


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
M2OldNmOptim<System,Agent,Model,ModelParam>::M2OldNmOptim(){
  name = "M2OldNm";
}

template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void M2OldNmOptim<System,Agent,Model,ModelParam>::
optim(System<Model,ModelParam> system,
      Agent<Model,ModelParam> & agent){

  qEval.preCompData(system.sD,system.fD);
  qEval.bellResFixData(system.sD,system.tD,system.fD,system.dD,
		       system.model,system.estParam);

  M2OldNmData<System,Agent,Model,ModelParam> qData;
  qData.qEval = qEval;
  qData.s = system;
  qData.a = agent;

  
  std::vector<double> par = agent.tp.getPar();
  int i,dim=agent.numFeatures;

  gsl_vector *x, *ss;
  x = gsl_vector_alloc(dim);
  for(i=0; i<dim; i++)
    gsl_vector_set(x,i,par.at(i));
  ss=gsl_vector_alloc(dim);
  gsl_vector_set_all(ss,1);

  gsl_multimin_function minex_func;
  minex_func.n=dim;
  minex_func.f=&M2OldNmObj<System,Agent,Model,ModelParam>;
  minex_func.params=&qData;

  const gsl_multimin_fminimizer_type *T=
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  s=gsl_multimin_fminimizer_alloc(T,dim);
  gsl_multimin_fminimizer_set(s,&minex_func,x,ss);

  double curSize;
  double size=.001;
  size_t iter=0;
  int status;
  
  do{
    iter++;
    status=gsl_multimin_fminimizer_iterate(s);
    if(status)
      break;
    curSize=gsl_multimin_fminimizer_size(s);
    status=gsl_multimin_test_size(curSize,size);

    // printf("iter % d: Q() = % 16.6f  ->  [",
    // 	   (int)iter,s->fval);
    // for(i=0; i<(dim-1); i++)
    //   printf(" % 10.6f,",gsl_vector_get(s->x,i));
    // printf(" % 10.6f ]\r",gsl_vector_get(s->x,i));
    // fflush(stdout);

  }while(status == GSL_CONTINUE && iter < 100);
  // std::cout << "\033[K";
  
  for(i=0; i<dim; i++)
    par.at(i) = gsl_vector_get(s->x,i);
  agent.tp.putPar(par);

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(ss);
}


template class M2OldNmEval<System,ProximalAgent,GravityModel,GravityParam>;

template class M2OldNmEval<System,RankAgent,GravityModel,GravityParam>;

template class M2OldNmEval<System,RankToyAgent,GravityModel,GravityParam>;


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
const int M2OldNmEval<System,Agent,Model,ModelParam>::numFeatures = 5;


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void M2OldNmEval<System,Agent,Model,ModelParam>::
preCompData(const SimData & sD, const FixedData & fD){
  int i,j;
  
  K=numFeatures + (numFeatures-1) + (numFeatures-1)*(numFeatures-2)/2;
  
  // obtain closest tp.numNeigh neighbors
  neighbors.clear();
  neighbors.resize(fD.numNodes);
  for(i=0; i<fD.numNodes; i++){
    std::priority_queue<std::pair<double,int> > pq;
    for(j=0; j<fD.numNodes; j++)
      if(i!=j)
	pq.push(std::pair<double,int>(-fD.dist.at(i*fD.numNodes+j),j));
    for(j=0; j<tp.numNeigh; j++){
      neighbors.at(i).push_back(pq.top().second);
      pq.pop();
    }
  }

  P = buildL2Pen(fD.numNodes);
  P += buildSpatialPen(sD,fD);
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void M2OldNmEval<System,Agent,Model,ModelParam>::
bellResFixData(const SimData & sD,
	       const TrtData & tD,
	       const FixedData & fD,
	       const DynamicData & dD,
	       const Model & m,
	       ModelParam & mP){
  int dim = (2*K-1)*fD.numNodes;
  // prep containers
  R.resize(dim);
  R.setZero();
  D0.resize(dim,dim);
  D0.setZero();
  Eigen::SparseMatrix<double> phi;
  Eigen::SparseMatrix<double> deltaQt;  
  Eigen::VectorXd Y(fD.numNodes);
  
  sD1T.clear();
  dD1T.clear();
  deltaQ.clear();
  
  SimData sDt;
  TrtData tDt;
  tDt.a.resize(fD.numNodes);
  tDt.p.resize(fD.numNodes);
  std::fill(tDt.a.begin(),tDt.a.end(),0);
  std::fill(tDt.p.begin(),tDt.p.end(),0);  

  std::vector<int> timeInft(fD.numNodes);
  std::fill(timeInft.begin(),timeInft.end(),0);

  std::vector< std::vector<int> > h;
  h = sD.history;

  std::vector<double> features;

  // setup initial SimData
  int i,t,I,status_i;
  for(t=0; t<sD.time; t++){

    // build complete history of simulation

    // clear containers and set simple things
    sDt.time=t;
    sDt.status = h.at(t);
    sDt.infected.clear();
    sDt.notInfec.clear();
    
    tDt.aPast=tDt.a;
    tDt.pPast=tDt.p;
    std::fill(tDt.a.begin(),tDt.a.end(),0);
    std::fill(tDt.p.begin(),tDt.p.end(),0);

    // set infected, notInfec, and treatments
    for(i=0; i<fD.numNodes; i++){
      status_i=sDt.status.at(i);
      if(status_i < 2){
	sDt.notInfec.push_back(i);
	if(status_i == 1)
	  tDt.p.at(i)=1;
      }
      else{
	sDt.infected.push_back(i);
	if(status_i == 3)
	  tDt.a.at(i)=1;
      }
    }

    // set timeInf
    sDt.numNotInfec = sDt.notInfec.size();
    sDt.numInfected = sDt.infected.size();

    for(i=0; i<sDt.numInfected; i++)
      timeInft.at(sDt.infected.at(i))++;
    sDt.timeInf=timeInft;

    // set newInf
    sDt.newInfec.clear();
    for(i=0; i<fD.numNodes; i++)
      if(sDt.timeInf.at(i)==1)
	sDt.newInfec.push_back(i);

    // set history
    sDt.history.clear();
    sDt.history.insert(sDt.history.end(),h.begin(),h.begin()+t);

    // add to containers
    if(t){ // t == 0 only needed for one time compute
      sD1T.push_back(sDt);
      dD1T.push_back(dD); // right now dD is completely empty ... trivial
    }


    // using the current values build D0
    features=getFeatures(sDt,tDt,fD,dD,m,mP);
    phi=featToPhi(features,fD.numNodes);
    
    deltaQt=phi.transpose();
    deltaQ.push_back(deltaQt);
    
    D0 += deltaQt * phi;

    Y.setZero();
    I=sDt.newInfec.size();
    for(i=0; i<I; i++)
      Y(sDt.newInfec.at(i))=-1;
    R += deltaQt * Y;
  }

  // sD and dD are at time T
  sD1T.push_back(sD);
  dD1T.push_back(dD);
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
Eigen::SparseMatrix<double> M2OldNmEval<System,Agent,Model,ModelParam>::
buildSpatialPen(const SimData & sD, const FixedData & fD){
  int i,j,k,node0,node1;
  Eigen::SparseMatrix<double> H,G,spatP;

  
  int dim=(2*K-1)*fD.numNodes;
  
  spatP.resize(dim,dim);
  spatP.setZero();
  
  // using the neighbor lookup, build the spatial penalty
  for(i=0; i<fD.numNodes; i++){
    node0=i;
    G.resize((2*K-1)*tp.numNeigh,dim);
    G.setZero();
    H.resize((2*K-1)*tp.numNeigh,dim);
    H.setZero();
    for(j=0; j<tp.numNeigh; j++){
      node1=neighbors.at(i).at(j);
      for(k=0; k<(2*K-1); k++){
	H.insert(j*(2*K-1)+k,node1*(2*K-1)+k)=1;
	G.insert(j*(2*K-1)+k,node0*(2*K-1)+k)=1;
      }
    }

    spatP += (H-G).transpose()*(H-G);
  }
  return spatP;
}



template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
Eigen::SparseMatrix<double> M2OldNmEval<System,Agent,Model,ModelParam>::
buildL2Pen(const int numNodes){
  Eigen::SparseMatrix<double> l2P((2*K-1)*numNodes,(2*K-1)*numNodes);
  int i,j;
  for(i=0; i<numNodes; i++)
    for(j=1; j<(2*K-1); j++)
      l2P.insert(i*(2*K-1)+j,i*(2*K-1)+j)=1;
  return l2P;
}





template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
std::vector<double> M2OldNmEval<System,Agent,Model,ModelParam>::
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

  // feature 0 (intercept)
  infFeat.col(featNum).ones();
  notFeat.col(featNum).ones();

  featNum++;
  

  // feature 1
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,0).t();
  
  featNum++;

  
  // feature 2
  SystemLight<Model,ModelParam> s(sD,tD,fD,dD,m,mP);
  std::vector<int> newInfec;
  int k,numNewInfec;
  for(i=0; i<tp.valReps; i++){
    s.reset();
    s.nextPoint();
    
    newInfec=s.sD.newInfec;
    s.nextPoint(1);
    
    newInfec.insert(newInfec.end(),s.sD.newInfec.begin(),s.sD.newInfec.end());
    
    numNewInfec = newInfec.size();
    for(j=0,itD0=newInfec.begin(); j<numNewInfec; j++,itD0++)
      for(k=0,itD1=sD.notInfec.begin(); k<sD.numNotInfec; k++,itD1++)
	if(*itD0 == *itD1)
	  notFeat(k,featNum)+=1.0/(double)tp.valReps;
  }

  infFeat.col(featNum) = (1.0-mP.infProbsSep) * notFeat.col(featNum);
  
  featNum++;


  // feature 3
  std::vector<double> infNoTrtLat,infNoTrtLong,notNoTrtLat,notNoTrtLong;
  int numInfNoTrt=0,numNotNoTrt=0;
  for(i=0; i<fD.numNodes; i++){
    if(sD.status.at(i) < 2 && tD.p.at(i) == 0){
      notNoTrtLat.push_back(fD.centroidsLat.at(i));
      notNoTrtLong.push_back(fD.centroidsLong.at(i));
      numNotNoTrt++;
    }
    else if(sD.status.at(i) >= 2 && tD.a.at(i) == 0){
      infNoTrtLat.push_back(fD.centroidsLat.at(i));
      infNoTrtLong.push_back(fD.centroidsLong.at(i));
      numInfNoTrt++;
    }
  }

  double minDist;
  for(i=0,itD0=sD.notInfec.begin(); i<sD.numNotInfec; i++,itD0++){
    minDist=0;
    for(j=0,itD1=sD.infected.begin(); j<sD.numInfected; j++,itD1++)
      if((1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1))) > minDist)
	minDist = 1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1));
    notFeat(i,featNum) = minDist/(halfPlaneDepth(fD.centroidsLong.at(*itD0),
						 fD.centroidsLat.at(*itD0),
						 numInfNoTrt,
						 infNoTrtLong,
						 infNoTrtLat)
				  + (1.0/(double)(numInfNoTrt+1)));
  }
  for(i=0,itD0=sD.infected.begin(); i<sD.numInfected; i++,itD0++){
    minDist=0;
    for(j=0,itD1=sD.notInfec.begin(); j<sD.numNotInfec; j++,itD1++)
      if((1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1))) > minDist)
	minDist = 1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1));
    infFeat(i,featNum) = minDist/(halfPlaneDepth(fD.centroidsLong.at(*itD0),
						 fD.centroidsLat.at(*itD0),
						 numNotNoTrt,
						 notNoTrtLong,
						 notNoTrtLat)
				  + (1.0/(double)(numNotNoTrt+1)));
  }
  featNum++;
  

  // feature 4
  std::vector<int>::const_iterator itD2,itD3;
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    notFeat(i,featNum) = std::log(1.0+totalDist);
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    infFeat(i,featNum) = std::log(1.0+totalDist);
  }

  featNum++;


  // put features into vector node by node
  std::vector<double> allFeat;
  arma::inplace_trans(notFeat);
  arma::inplace_trans(infFeat);
  for(i=0,j=0,k=0; i<fD.numNodes; i++){
    if(sD.status.at(i) < 2){
      allFeat.insert(allFeat.end(),notFeat.begin()+j*numFeatures,
		     notFeat.begin()+(j+1)*numFeatures);
      j++;
    }
    else{
      allFeat.insert(allFeat.end(),infFeat.begin()+k*numFeatures,
		     infFeat.begin()+(k+1)*numFeatures);
      k++;
    }
  }
  
  return allFeat;
}





template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
inline Eigen::SparseMatrix<double> M2OldNmEval<System,Agent,Model,ModelParam>::
featToPhi(const std::vector<double> & feat, const int numNodes){

  int n,i,j,dim=(2*K-1)*numNodes;

  // create full features with interactions
  std::vector<double> featFull;
  featFull.reserve(K*numNodes);
  for(n=0; n<numNodes; n++){
    for(i=0; i<numFeatures; i++){
      featFull.push_back(feat.at(n*numFeatures +i));
    }

    for(i=1; i<numFeatures; i++){
      for(j=i; j<numFeatures; j++){
	featFull.push_back(feat.at(n*numFeatures + i)*
			   feat.at(n*numFeatures + j));
      }
    }
  }

  
  // enter in each nodes features
  Eigen::SparseMatrix<double> Di(numNodes,dim);
  for(n=0; n<numNodes; n++)
    for(i=0; i<K; i++)
      Di.insert(n,(2*K-1)*n + i) = featFull.at(n*K + i);
	 

  // enter in average of neighbors features
  double neighAvg;
  int node1;
  for(n=0; n<numNodes; n++){
    for(i=0; i<(K-1); i++){
      neighAvg=0.0;
      for(j=0; j<tp.numNeigh; j++){
	node1=neighbors.at(n).at(j);
	neighAvg+=featFull.at(node1*K + 1 + i);
      }
      neighAvg/=(double)tp.numNeigh;
      Di.insert(n,(2*K-1)*n + K + i) = neighAvg;
    }
  }
    
  return Di;
}





template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void M2OldNmEval<System,Agent,Model,ModelParam>::
bellResPolData(const int time,
	       const FixedData & fD,
	       const Model & m,
	       ModelParam & mP,
	       Agent<Model,ModelParam> a){

  int dim = (2*K-1)*fD.numNodes;
  Eigen::SparseMatrix<double> phi;

  D1.resize(dim,dim);
  D1.setZero();

  std::vector<double> features;
  TrtData tDt;
  tDt.a.resize(fD.numNodes);
  tDt.p.resize(fD.numNodes);
  tDt.aPast.resize(fD.numNodes);
  tDt.pPast.resize(fD.numNodes);
  std::fill(tDt.a.begin(),tDt.a.end(),0);
  std::fill(tDt.p.begin(),tDt.p.end(),0);  
  std::fill(tDt.aPast.begin(),tDt.aPast.end(),0);
  std::fill(tDt.pPast.begin(),tDt.pPast.end(),0);  
  
  int t,i,j;
  for(t=0; t<time; t++){
    std::fill(tDt.a.begin(),tDt.a.end(),0);
    std::fill(tDt.p.begin(),tDt.p.end(),0);  
    std::fill(tDt.aPast.begin(),tDt.aPast.end(),0);
    std::fill(tDt.pPast.begin(),tDt.pPast.end(),0);
    
    for(i=0; i<fD.numNodes; i++){
      if(sD1T.at(t).history.at(t).at(i)==1)
	tDt.pPast.at(i)=1;
      else if(sD1T.at(t).history.at(t).at(i)==3)
	tDt.aPast.at(i)=1;
    }

    if((t+1)>=fD.trtStart){
      phi.resize(fD.numNodes,fD.numNodes*(2*K-1));
      phi.setZero();      
      for(j=0; j<tp.polReps; j++){
	std::fill(tDt.a.begin(),tDt.a.end(),0);
	std::fill(tDt.p.begin(),tDt.p.end(),0);  
	a.applyTrt(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
    
	features = getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	phi += featToPhi(features,fD.numNodes);
      }
      phi /= (double)tp.polReps;
    }
    else{
      features = getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      phi = featToPhi(features,fD.numNodes);
    }

    D1 += deltaQ.at(t) * phi;
  }

  D1 *= tp.gamma; // discount factor

  D = D1 - D0;
}

template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
void M2OldNmEval<System,Agent,Model,ModelParam>::
solve(){
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > simplicialLDLT;
  simplicialLDLT.compute(D.transpose()*D + tp.lambda*P);

  beta = simplicialLDLT.solve(-D.transpose()*R);
}

template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double M2OldNmEval<System,Agent,Model,ModelParam>::
bellRes(){
  return (R + D*beta).squaredNorm();
}


template < template<typename,typename> class System,
	   template<typename,typename> class Agent,
	   class Model,class ModelParam>
double M2OldNmEval<System,Agent,Model,ModelParam>::
qFn(const SimData & sD,
    TrtData & tD,
    const FixedData & fD,
    const DynamicData & dD,
    const Model & m,
    ModelParam & mP,
    Agent<Model,ModelParam> a){
  // now evaluate Q-function
  Eigen::SparseMatrix<double> phi;
  phi.resize(fD.numNodes,fD.numNodes*(2*K-1));
  phi.setZero();
  int j;
  std::vector<double> features;
  for(j=0; j<tp.polReps; j++){
    std::fill(tD.a.begin(),tD.a.end(),0);
    std::fill(tD.p.begin(),tD.p.end(),0);
    a.applyTrt(sD,tD,fD,dD,m,mP);
    
    features = getFeatures(sD,tD,fD,dD,m,mP);
    phi += featToPhi(features,fD.numNodes);
  }
  phi /= (double)tp.polReps;

  return (phi*beta).sum();
}
