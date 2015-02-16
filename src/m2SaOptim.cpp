#include "m2SaOptim.hpp"



std::vector<double> M2SaEvalTunePar::getPar() const{
  return std::vector<double> (0);
}



void M2SaEvalTunePar::putPar(const std::vector<double> & par){
}



template class M2SaOptim<System<GravityModel,GravityParam,
				GravityModel,GravityParam>,
			 RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
				      GravityModel,GravityParam>,
			 FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
				     GravityModel,GravityParam>,
			 GravityModel,GravityParam>;

template class M2SaOptim<System<RangeModel,RangeParam,
				RangeModel,RangeParam>,
			 RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				      RangeModel,RangeParam>,
			 FeaturesInt<ToyFeatures2<RangeModel,RangeParam>,
				     RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class M2SaOptim<System<CaveModel,CaveParam,
				CaveModel,CaveParam>,
			 RankToyAgent<ToyFeatures2<CaveModel,CaveParam>,
				      CaveModel,CaveParam>,
			 FeaturesInt<ToyFeatures2<CaveModel,CaveParam>,
				     CaveModel,CaveParam>,
			 CaveModel,CaveParam>;

template class M2SaOptim<System<GravityModel,GravityParam,
				RangeModel,RangeParam>,
			 RankToyAgent<ToyFeatures2<RangeModel,RangeParam>,
				      RangeModel,RangeParam>,
			 FeaturesInt<ToyFeatures2<RangeModel,RangeParam>,
				     RangeModel,RangeParam>,
			 RangeModel,RangeParam>;

template class M2SaOptim<System<GravityModel,GravityParam,
				CaveModel,CaveParam>,
			 RankToyAgent<ToyFeatures2<CaveModel,CaveParam>,
				      CaveModel,CaveParam>,
			 FeaturesInt<ToyFeatures2<CaveModel,CaveParam>,
				     CaveModel,CaveParam>,
			 CaveModel,CaveParam>;



template <class S, class A, class F,
	  class M,class MP>
M2SaOptim<S,A,F,M,MP>::M2SaOptim(){
  name = "M2Sa";
}



template <class S, class A, class F,
	  class M,class MP>
void M2SaOptim<S,A,F,M,MP>::
optim(const S & system,
      A & agent){

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);

  qEval.preCompData(s.sD,s.fD);
  qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,
		       s.modelEst,s.paramEst);
  


  double sd = qEval.tp.sdStart;
  std::vector<double> w;
  std::vector<double> bestW;
  double q,bestQ;
  
  qEval.bellResPolData(s.sD.time,s.fD,
		       s.modelEst,s.paramEst,agent);
  qEval.solve();
  
  bestW = agent.tp.getPar();
  bestQ=qEval.qFn(s.sD,s.tD,s.fD,s.dD,
		  s.modelEst,s.paramEst,agent);
  
  int i;  
  while(sd > qEval.tp.sdStop){
    w.clear();
    for(i=0; i<agent.f.numFeatures; i++)
      w.push_back(bestW.at(i) +sd*njm::rnorm01());
    njm::l2norm(w);
    agent.tp.putPar(w);

    qEval.bellResPolData(s.sD.time,s.fD,
			 s.modelEst,s.paramEst,agent);
    qEval.solve();
    q=qEval.qFn(s.sD,s.tD,s.fD,s.dD,
		s.modelEst,s.paramEst,agent);
    
    if(q > bestQ){
      bestQ = q;
      bestW = w;
      sd*=qEval.tp.sdJump;
    }
    else{
      sd*=qEval.tp.sdDecay;
    }
  }

  agent.tp.putPar(bestW);
}


template class M2SaEval<System<GravityModel,GravityParam,
			       GravityModel,GravityParam>,
			RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
				     GravityModel,GravityParam>,
			FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
				    GravityModel,GravityParam>,
			GravityModel,GravityParam>;

template <class S, class A, class F,
	  class M,class MP>
M2SaEval<S,A,F,M,MP>::M2SaEval(){
  tp.polReps = 10;
  // radius = 50;
  tp.numNeigh = 5;

  tp.gamma = .95;
  tp.lambda = 3000.0;
  
  tp.jitter = .1;
  tp.tol = .005;
  tp.rate = 5;
  tp.rateDecay = .975;

  tp.sdStart = 2.0;
  tp.sdStop = 0.5;
  tp.sdJump = 5.0;
  tp.sdDecay = 0.975;
}



template <class S, class A, class F,
	  class M,class MP>
void M2SaEval<S,A,F,M,MP>::
preCompData(const SimData & sD, const FixedData & fD){
  int i,j;
  
  K=f.numFeatures + (f.numFeatures-1) + (f.numFeatures-1)*(f.numFeatures-2)/2;
  
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



template <class S, class A, class F,
	  class M,class MP>
void M2SaEval<S,A,F,M,MP>::
bellResFixData(const SimData & sD,
	       const TrtData & tD,
	       const FixedData & fD,
	       const DynamicData & dD,
	       const M & m,
	       MP & mP){
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
    f.preCompData(sDt,tDt,fD,dD,m,mP);
    f.getFeatures(sDt,tDt,fD,dD,m,mP);
    features=feat2Vec(fD.numNodes,sDt.status);
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



template <class S, class A, class F,
	  class M,class MP>
Eigen::SparseMatrix<double> M2SaEval<S,A,F,M,MP>::
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



template <class S, class A, class F,
	  class M,class MP>
Eigen::SparseMatrix<double> M2SaEval<S,A,F,M,MP>::
buildL2Pen(const int numNodes){
  Eigen::SparseMatrix<double> l2P((2*K-1)*numNodes,(2*K-1)*numNodes);
  int i,j;
  for(i=0; i<numNodes; i++)
    for(j=1; j<(2*K-1); j++)
      l2P.insert(i*(2*K-1)+j,i*(2*K-1)+j)=1;
  return l2P;
}




template <class S, class A, class F,
	  class M,class MP>
std::vector<double> M2SaEval<S,A,F,M,MP>::
feat2Vec(const int numNodes,
	 const std::vector<int> & status){

  // put features into vector node by node
  std::vector<double> allFeat;
  arma::inplace_trans(f.notFeat);
  arma::inplace_trans(f.infFeat);
  int i,j,k;
  for(i=0,j=0,k=0; i<numNodes; i++){
    if(status.at(i) < 2){
      allFeat.insert(allFeat.end(),f.notFeat.begin()+j*f.numFeatures,
		     f.notFeat.begin()+(j+1)*f.numFeatures);
      j++;
    }
    else{
      allFeat.insert(allFeat.end(),f.infFeat.begin()+k*f.numFeatures,
		     f.infFeat.begin()+(k+1)*f.numFeatures);
      k++;
    }
  }

  arma::inplace_trans(f.notFeat);
  arma::inplace_trans(f.infFeat);
  
  return allFeat;
}





template <class S, class A, class F,
	  class M,class MP>
inline Eigen::SparseMatrix<double>
M2SaEval<S,A,F,M,MP>::
featToPhi(const std::vector<double> & feat, const int numNodes){

  int n,i,j,dim=(2*K-1)*numNodes;

  // create full features with interactions
  std::vector<double> featFull;
  featFull.reserve(K*numNodes);
  for(n=0; n<numNodes; n++){
    for(i=0; i<f.numFeatures; i++){
      featFull.push_back(feat.at(n*f.numFeatures +i));
    }

    for(i=1; i<f.numFeatures; i++){
      for(j=i; j<f.numFeatures; j++){
	featFull.push_back(feat.at(n*f.numFeatures + i)*
			   feat.at(n*f.numFeatures + j));
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





template <class S, class A, class F,
	  class M,class MP>
void M2SaEval<S,A,F,M,MP>::
bellResPolData(const int time,
	       const FixedData & fD,
	       const M & m,
	       MP & mP,
	       A a){

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

	f.preCompData(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	f.getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	features = feat2Vec(fD.numNodes,sD1T.at(t).status);
	phi += featToPhi(features,fD.numNodes);
      }
      phi /= (double)tp.polReps;
    }
    else{
      f.preCompData(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      f.getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      features = feat2Vec(fD.numNodes,sD1T.at(t).status);
      phi = featToPhi(features,fD.numNodes);
    }

    D1 += deltaQ.at(t) * phi;
  }

  D1 *= tp.gamma; // discount factor

  D = D1 - D0;
}



template <class S, class A, class F,
	  class M,class MP>
void M2SaEval<S,A,F,M,MP>::
solve(){
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > simplicialLDLT;
  simplicialLDLT.compute(D.transpose()*D + tp.lambda*P);

  beta = simplicialLDLT.solve(-D.transpose()*R);
}



template <class S, class A, class F,
	  class M,class MP>
double M2SaEval<S,A,F,M,MP>::
bellRes(){
  return (R + D*beta).squaredNorm();
}



template <class S, class A, class F,
	  class M,class MP>
double M2SaEval<S,A,F,M,MP>::
qFn(const SimData & sD,
    TrtData & tD,
    const FixedData & fD,
    const DynamicData & dD,
    const M & m,
    MP & mP,
    A a){
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

    f.preCompData(sD,tD,fD,dD,m,mP);
    f.getFeatures(sD,tD,fD,dD,m,mP);
    features = feat2Vec(fD.numNodes,sD.status);
    phi += featToPhi(features,fD.numNodes);
  }
  phi /= (double)tp.polReps;

  return (phi*beta).sum();
}
