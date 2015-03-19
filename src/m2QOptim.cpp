#include "m2QOptim.hpp"



std::vector<double> M2QEvalTunePar::getPar() const{
  return std::vector<double> (0);
}



void M2QEvalTunePar::putPar(const std::vector<double> & par){
}



template class M2QOptim<System<GravityModel,GravityParam,
				GravityModel,GravityParam>,
			 RankAgent<ToyFeatures2<GravityModel,GravityParam>,
				      GravityModel,GravityParam>,
			 FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
				     GravityModel,GravityParam>,
			 GravityModel,GravityParam>;


template class
M2QOptim<System<GravityTimeInfExpCavesModel,
		GravityTimeInfExpCavesParam,
		GravityTimeInfExpCavesModel,
		GravityTimeInfExpCavesParam>,
	 RankAgent<ToyFeatures2<GravityTimeInfExpCavesModel,
				GravityTimeInfExpCavesParam>,
		   GravityTimeInfExpCavesModel,
		   GravityTimeInfExpCavesParam>,
	 FeaturesInt<ToyFeatures2<GravityTimeInfExpCavesModel,
				  GravityTimeInfExpCavesParam>,
		     GravityTimeInfExpCavesModel,GravityTimeInfExpCavesParam>,
	 GravityTimeInfExpCavesModel,GravityTimeInfExpCavesParam>;



template <class S, class A, class F,
	  class M,class MP>
M2QOptim<S,A,F,M,MP>::M2QOptim(){
  name = "M2Q";
}



template <class S, class A, class F,
	  class M,class MP>
void M2QOptim<S,A,F,M,MP>::
optim(const S & system,
      A & agent){

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);

  qEval.preCompData(s.sD,s.fD);
  qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,
		       s.modelEst,s.paramEst);
  


  // double sd = qEval.tp.sdStart;
  // std::vector<double> w;
  // std::vector<double> bestW;
  // double q,bestQ;
  
  // qEval.bellResPolData(s.sD.time,s.fD,
  // 		       s.modelEst,s.paramEst,agent);
  // qEval.solve();
  
  // bestW = agent.tp.getPar();
  // bestQ=qEval.qFn(s.sD,s.tD,s.fD,s.dD,
  // 		  s.modelEst,s.paramEst,agent);
  
  // int i;  
  // while(sd > qEval.tp.sdStop){
  //   w.clear();
  //   for(i=0; i<agent.f.numFeatures; i++)
  //     w.push_back(bestW.at(i) +sd*njm::rnorm01());
  //   njm::l2norm(w);
  //   agent.tp.putPar(w);

  //   qEval.bellResPolData(s.sD.time,s.fD,
  // 			 s.modelEst,s.paramEst,agent);
  //   qEval.solve();
  //   q=qEval.qFn(s.sD,s.tD,s.fD,s.dD,
  // 		s.modelEst,s.paramEst,agent);
    
  //   if(q > bestQ){
  //     bestQ = q;
  //     bestW = w;
  //     sd*=qEval.tp.sdJump;
  //   }
  //   else{
  //     sd*=qEval.tp.sdDecay;
  //   }
  // }

  // agent.tp.putPar(bestW);
}


template class M2QEval<System<GravityModel,GravityParam,
			       GravityModel,GravityParam>,
			RankAgent<ToyFeatures2<GravityModel,GravityParam>,
				     GravityModel,GravityParam>,
			FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
				    GravityModel,GravityParam>,
			GravityModel,GravityParam>;



template class M2QEval<System<GravityTimeInfExpCavesModel,
			      GravityTimeInfExpCavesParam,
			       GravityTimeInfExpCavesModel,
			      GravityTimeInfExpCavesParam>,
			RankAgent<ToyFeatures2<GravityTimeInfExpCavesModel,
					       GravityTimeInfExpCavesParam>,
				     GravityTimeInfExpCavesModel,
				  GravityTimeInfExpCavesParam>,
			FeaturesInt<ToyFeatures2<GravityTimeInfExpCavesModel,
						 GravityTimeInfExpCavesParam>,
				    GravityTimeInfExpCavesModel,
				    GravityTimeInfExpCavesParam>,
			GravityTimeInfExpCavesModel,
		       GravityTimeInfExpCavesParam>;




template <class S, class A, class F,
	  class M,class MP>
M2QEval<S,A,F,M,MP>::M2QEval(){
  tp.polReps = 10;
  tp.numNeigh = 5;

  tp.gamma = .95;
  // tp.lambda = 3000.0;
  tp.lambda = 1000;

  tp.dfLat = 10;
  tp.dfLong = 10;

  tp.bootReps = 10;
  tp.bootSize = 0.6;
}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
preCompData(const SimData & sD, const FixedData & fD){
  int i,j;

  numNodes = fD.numNodes;
  
  // explanation is in header
  numFeat = f.numFeatures + (f.numFeatures-1) +
    (f.numFeatures-1)*(f.numFeatures-2)/2;
  lenPsi = numFeat*2 - 1;
  
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


  // compute the bspline matrix with MDS coordinates
  const int deg = 4; // (k - 1) polynomial spline
  const int nbreakLat = tp.dfLat + 2 - deg;
  const int nbreakLong = tp.dfLong + 2 - deg;


  gsl_bspline_workspace *bsLat, *bsLong;
  bsLat = gsl_bspline_alloc(deg,nbreakLat);
  bsLong = gsl_bspline_alloc(deg,nbreakLong);


  double latMin = *std::min_element(fD.centroidsMdsLat.begin(),
				    fD.centroidsMdsLat.end());
  double latMax = *std::max_element(fD.centroidsMdsLat.begin(),
				    fD.centroidsMdsLat.end());
  double longMin = *std::min_element(fD.centroidsMdsLong.begin(),
				     fD.centroidsMdsLong.end());
  double longMax = *std::max_element(fD.centroidsMdsLong.begin(),
				     fD.centroidsMdsLong.end());

  
  gsl_bspline_knots_uniform(latMin,latMax,bsLat);
  gsl_bspline_knots_uniform(longMin,longMax,bsLong);

  gsl_vector *bLat, *bLong;
  bLat = gsl_vector_alloc(tp.dfLat);
  bLong = gsl_vector_alloc(tp.dfLong);


  phiL.clear();
  int u,v;
  double bb;
  for(i = 0; i < fD.numNodes; ++i){
    Eigen::SparseMatrix<double> phi(lenPsi*tp.dfLat*tp.dfLong,lenPsi);
    phi.reserve(lenPsi*tp.dfLat*tp.dfLong);
    
    gsl_bspline_eval(fD.centroidsMdsLat.at(i),
		     bLat,bsLat);
    gsl_bspline_eval(fD.centroidsMdsLong.at(i),
		     bLong,bsLong);
    for(j = 0; j < lenPsi; ++j){
      for(u = 0; u < tp.dfLat; ++u){
	for(v = 0; v < tp.dfLong; ++v){
	  bb = gsl_vector_get(bLat,u)*gsl_vector_get(bLong,v);
	  if(bb != 0.0)
	    phi.insert(j*tp.dfLat*tp.dfLong + u*tp.dfLong + v,j) = bb;
	}
      }
    }
    phiL.push_back(phi);
  }

  gsl_bspline_free(bsLat);
  gsl_bspline_free(bsLong);
  gsl_vector_free(bLat);
  gsl_vector_free(bLong);
}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
bellResFixData(const SimData & sD,
	       const TrtData & tD,
	       const FixedData & fD,
	       const DynamicData & dD,
	       const M & m,
	       MP & mP){
  int dim = tp.dfLat*tp.dfLong*lenPsi;
  // prep containers
  RL.resize(fD.numNodes);
  std::fill(RL.begin(),RL.end(),Eigen::VectorXd(dim));

  D0L.resize(fD.numNodes);
  std::fill(D0L.begin(),D0L.end(),Eigen::SparseMatrix<double>(dim,dim));

  // make sure containers are zero'd out
  int i;
  for(i = 0; i < fD.numNodes; ++i){
    RL.at(i).setZero();
    D0L.at(i).setZero();
  }
  
  sD1T.clear();
  dD1T.clear();
  
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

  std::vector<Eigen::SparseMatrix<double> > phiPsiL;

  // setup initial SimData
  int t,status_i,numNewInfec;
  for(t=0; t<sD.time; t++){

    njm::timer.start("fixBuild");
    
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
    if(t > 0){ // t == 0 only needed for one time compute
      sD1T.push_back(sDt);
      dD1T.push_back(dD); // right now dD is completely empty ... trivial
    }

    njm::timer.stop("fixBuild");    

    njm::timer.start("fixFeat");    
    // using the current values build D0
    f.preCompData(sDt,tDt,fD,dD,m,mP);
    f.getFeatures(sDt,tDt,fD,dD,m,mP);
    features=feat2Vec(fD.numNodes,sDt.status);

    njm::timer.stop("fixFeat");    

    njm::timer.start("fixPhiPsi");
    phiPsiL=featToPhiPsi(features,fD.numNodes);
    
    phiPsiTL.push_back(phiPsiL);
    njm::timer.stop("fixPhiPsi");

    for(i = 0; i < fD.numNodes; ++i){
      njm::timer.start("fixD0");
      D0L.at(i) += phiPsiL.at(i) * phiPsiL.at(i).transpose();
      njm::timer.stop("fixD0");
    }

    njm::timer.start("fixR");
    numNewInfec = sDt.newInfec.size();
    for(i = 0; i < numNewInfec; ++i)
      RL.at(i) += phiPsiL.at(sDt.newInfec.at(i)) * (1.0/double(fD.numNodes));
    njm::timer.stop("fixR");
  }

  // sD and dD are at time T
  sD1T.push_back(sD);
  dD1T.push_back(dD);
}




template <class S, class A, class F,
	  class M,class MP>
std::vector<double> M2QEval<S,A,F,M,MP>::
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
inline std::vector<Eigen::SparseMatrix<double> >
M2QEval<S,A,F,M,MP>::
featToPhiPsi(const std::vector<double> & feat, const int numNodes){

  int n,i,j;

  // create full features with interactions
  std::vector<double> featFull;
  featFull.reserve(numFeat*numNodes);
  for(n=0; n<numNodes; n++){
    for(i=0; i<f.numFeatures; i++){
      featFull.push_back(feat.at(n*f.numFeatures +i));
    }

    for(i=1; i<f.numFeatures; i++){// start at 1 since 0 is the intercept
      for(j=i; j<f.numFeatures; j++){
	featFull.push_back(feat.at(n*f.numFeatures + i)*
			   feat.at(n*f.numFeatures + j));
      }
    }
  }

  // for each location create psi vector
  double avg;
  std::vector<Eigen::SparseMatrix<double> > mats;
  for(n = 0; n < numNodes; ++n){
    Eigen::SparseMatrix<double> psiL(lenPsi,1);
    
    for(i = 0; i < numFeat; ++i)
      psiL.insert(i,0) = featFull.at(n*numFeat + i);
    
    for(i = 1; i < numFeat; ++i){
      avg = 0;
      for(j = 0; j < tp.numNeigh; ++j)
	avg += featFull.at(neighbors.at(n).at(j)*numFeat + i);
      avg /= double(tp.numNeigh);
      psiL.insert(numFeat + i - 1,0) = avg;
    }
    
    mats.push_back(phiL.at(n) * psiL);
  }
  
  return mats;
}





template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
bellResPolData(const int time,
	       const FixedData & fD,
	       const M & m,
	       MP & mP,
	       A a){

  int dim = tp.dfLat*tp.dfLong*lenPsi;
  std::vector<Eigen::SparseMatrix<double> > phiPsiL;
  std::vector<Eigen::SparseMatrix<double> > phiPsiLavg;

  D1L.resize(fD.numNodes);
  std::fill(D1L.begin(),D1L.end(),Eigen::SparseMatrix<double>(dim,dim));

  int i;
  for(i = 0; i < fD.numNodes; ++i)
    D1L.at(i).setZero();
  

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
  
  int t,j,k;
  for(t=0; t<time; t++){
    njm::timer.start("polBuild");
    
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

    njm::timer.stop("polBuild");

    if((t+1)>=fD.trtStart){
      for(j=0; j<tp.polReps; j++){
	njm::timer.start("polFeat");
	std::fill(tDt.a.begin(),tDt.a.end(),0);
	std::fill(tDt.p.begin(),tDt.p.end(),0);  
	a.applyTrt(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);

	f.preCompData(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	f.getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	features = feat2Vec(fD.numNodes,sD1T.at(t).status);

	njm::timer.stop("polFeat");

	njm::timer.start("polPhiPsi");
	phiPsiL = featToPhiPsi(features,fD.numNodes);
	if(j == 0)
	  phiPsiLavg = phiPsiL;
	else
	  for(k = 0; k < fD.numNodes; ++k)
	    phiPsiLavg.at(k) += phiPsiL.at(k);
	njm::timer.stop("polPhiPsi");
      }
      njm::timer.start("polPhiPsi");
      for(k = 0; k < fD.numNodes; ++k)
	phiPsiLavg.at(k) /= double(tp.polReps);
      njm::timer.stop("polPhiPsi");
    }
    else{
      njm::timer.start("polFeat");
      f.preCompData(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      f.getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      features = feat2Vec(fD.numNodes,sD1T.at(t).status);
      njm::timer.stop("polFeat");

      njm::timer.start("polPhiPsi");
      phiPsiLavg = featToPhiPsi(features,fD.numNodes);
      njm::timer.stop("polPhiPsi");
    }

    njm::timer.start("polD1");
    for(k = 0; k < fD.numNodes; ++k)
      D1L.at(k) += phiPsiTL.at(t).at(k) * phiPsiLavg.at(k).transpose();
    njm::timer.stop("polD1");
  }

  njm::timer.start("polFinish");

  for(k = 0; k < fD.numNodes; ++k)
    D1L.at(k) *= tp.gamma; // discount factor

  njm::timer.stop("polFinish");
}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
solve(){
  njm::timer.start("solve");
  
  Eigen::SparseMatrix<double> P(tp.dfLat*tp.dfLong*lenPsi,
				tp.dfLat*tp.dfLong*lenPsi);
  P.setIdentity();

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  solver.compute(D.transpose() * D + tp.lambda*P);

  if(solver.info() != Eigen::Success){
    std::cout << "In M2QEval::solve(): decomposition failed."
	      << std::endl;
    throw(1);
  }

  beta = solver.solve(-D.transpose() * R);

  njm::timer.start("stop");  
}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
buildRD(){
  std::vector<int> nodes;
  nodes.resize(numNodes);
  int i = 0;
  std::for_each(nodes.begin(),nodes.end(),
		[&i](int & x){
		  x = i++;
		});
  
  buildRD(nodes);
}


template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
buildRD(const std::vector<int> nodes){
  njm::timer.start("buildRD");
  
  int i,I = nodes.size();
  int dim = tp.dfLat*tp.dfLong*lenPsi;  
  R.resize(dim);
  D0.resize(dim,dim);
  D1.resize(dim,dim);

  R.setZero();
  D0.setZero();
  D1.setZero();
  for(i = 0; i < I; ++i)
    D0 += D0L.at(i);
  for(i = 0; i < I; ++i)
    D1 += D1L.at(i);
  for(i = 0; i < I; ++i)
    R += RL.at(i);

  D = D1 - D0;

  njm::timer.stop("buildRD");
}





template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
buildRD1(){
  std::vector<int> nodes;
  nodes.resize(numNodes);
  int i = 0;
  std::for_each(nodes.begin(),nodes.end(),
		[&i](int & x){
		  x = i++;
		});
  
  buildRD1(nodes);
}


template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
buildRD1(const std::vector<int> nodes){
  njm::timer.start("buildRD1");
  
  int i,I = nodes.size();
  int dim = tp.dfLat*tp.dfLong*lenPsi;  
  D1.resize(dim,dim);

  D1.setZero();
  for(i = 0; i < I; ++i)
    R += RL.at(i);

  D = D1 - D0;

  njm::timer.stop("buildRD1");
}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
tune(){
  int i,j,b,bS = (tp.bootSize * numNodes + 1);
  std::vector<int> nodes;
  ind.reserve(numNodes);
  for(i = 0; i < numNodes; ++i)
    nodes.push_back(i);

  int numLambdaPows = 10;
  std::vector<double> lambdaPows;
  for(i = 0; i < numLambdaPows; ++i)
    lambdaPows.push_back(i+1);

  std::vector<double> lambdaCV;
  lambdaCV.resize(numLambdaVals);
  std::fill(lambdaCV.begin(),lambda.end(),0.0);


  
  std::pair<double,int> top;
  for(b = 0; b < tp.bootReps; ++b){
    std::priority_queue<std::pair<double,int> > ordNodes;
    for(i = 0; i < numNodes; ++i)
      ordNodes.push(std::pair<double,int>(njm::runif01(),nodes.at(i)));


    std::vector<int> selNodes;
    for(i = 0; i < bS; ++i){
      top = ordNodes.top();
      ordNodes.pop();
      selNodes.push_back(top.second());
    }

    buildRD(selNodes);

    for(i = 0; i < lambdaVals; ++i){
      tp.lambda = std::pow(2.0,lambdaPows.at(i));
      
      solve();

      lambdaCV.at(i) += bellRes();
    }
  }


  // pick off the best
  std::priority_queue<std::pair<double,double> > ordLambda;
  for(i = 0; i < lambdaPows; ++i)
    ordLambda_push(std::pair<double,double>(-lambdaCV.at(i),lambdaPows.at(i)));
  double bestPow = ordLambda.top().second;

  // now grid it up finer
  lambdaPows.clear();
  lambdaCV.clear();
  
  numLambdaPows = 10 + 1;
  for(i = 0; i < numLambdaPows; ++i)
    if(i != 5)
      lambdaPows.push_back(bestPow - double(i - 5)*0.2);
  --numLambdaPows;


  // error for new lambdas
  for(b = 0; b < tp.bootReps; ++b){
    std::priority_queue<std::pair<double,int> > ordNodes;
    for(i = 0; i < numNodes; ++i)
      ordNodes.push(std::pair<double,int>(njm::runif01(),nodes.at(i)));


    std::vector<int> selNodes;
    for(i = 0; i < bS; ++i){
      top = ordNodes.top();
      ordNodes.pop();
      selNodes.push_back(top.second());
    }

    buildRD(selNodes);

    for(i = 0; i < lambdaVals; ++i){
      tp.lambda = std::pow(2.0,lambdaPows.at(i));
      
      solve();

      lambdaCV.at(i) += bellRes();
    }
  }


  for(i = 0; i < lambdaPows; ++i)
    ordLambda_push(std::pair<double,double>(-lambdaCV.at(i),lambdaPows.at(i)));
  double bestPow = ordLambda.top().second;
  
}




template <class S, class A, class F,
	  class M,class MP>
double M2QEval<S,A,F,M,MP>::
bellRes(){
  return (R + D*beta).squaredNorm();
}



template <class S, class A, class F,
	  class M,class MP>
double M2QEval<S,A,F,M,MP>::
qFn(const SimData & sD,
    TrtData & tD,
    const FixedData & fD,
    const DynamicData & dD,
    const M & m,
    MP & mP,
    A a){

  njm::timer.start("qFn");
  
  // now evaluate Q-function
  std::vector<Eigen::SparseMatrix<double> > phiPsiL;
  Eigen::SparseMatrix<double> phiPsi;
  phiPsi.resize(tp.dfLat*tp.dfLong*lenPsi,1);
  phiPsi.setZero();
  int j,k;
  std::vector<double> features;
  for(j=0; j<tp.polReps; j++){
    std::fill(tD.a.begin(),tD.a.end(),0);
    std::fill(tD.p.begin(),tD.p.end(),0);
    a.applyTrt(sD,tD,fD,dD,m,mP);

    f.preCompData(sD,tD,fD,dD,m,mP);
    f.getFeatures(sD,tD,fD,dD,m,mP);
    features = feat2Vec(fD.numNodes,sD.status);
    phiPsiL = featToPhiPsi(features,fD.numNodes);

    for(k = 0; k < fD.numNodes; ++k)
      phiPsi += phiPsiL.at(k);
  }
  phiPsi /= (double)tp.polReps;

  njm::timer.stop("qFn");
  
  return (phiPsi.transpose()*beta).sum();
}
