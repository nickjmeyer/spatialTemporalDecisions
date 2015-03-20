#include "m2QOptim.hpp"



std::vector<double> M2QEvalTunePar::getPar() const{
  return std::vector<double> (0);
}



void M2QEvalTunePar::putPar(const std::vector<double> & par){
}


std::vector<double> M2QOptimTunePar::getPar() const{
  return std::vector<double> (0);
}



void M2QOptimTunePar::putPar(const std::vector<double> & par){
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
	  class M, class MP>
M2QOptim<S,A,F,M,MP>::M2QOptim(){
  name = "M2Q";


  tp.C = 5.0;

  tp.t = 1.0;

  tp.ell = 1.75;

  tp.muMin = 0.1;

  tp.A = 50;
  tp.B = 1;
}


template <class S, class A, class F,
	  class M, class MP>
void M2QOptim<S,A,F,M,MP>::reset(){
}


template <class S, class A, class F,
	  class M,class MP>
void M2QOptim<S,A,F,M,MP>::
optim(const S & system,
      A & agent){

  // first three steps use weights {1,1,...}
  if(system.sD.time < (system.fD.trtStart + 3))
    return;

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);  

  qEval.preCompData(s.sD,s.fD);
  qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,
		       s.modelEst,s.paramEst);
  qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,agent);

  if(system.sD.time == (system.fD.trtStart + 3))
    qEval.tune(system.sD.status);
  
  qEval.buildRD();
  


  std::vector<double> par=agent.tp.getPar();
  int i,converged=0,numPar = par.size();
  std::vector<double> h(numPar,0.0);
  std::vector<double> parPH(numPar,0.0),parMH(numPar,0.0);


  double valP, valM;
  
  int iter=1;
  
  double mu = tp.A/std::pow(tp.B+iter,tp.ell);
  double cm = tp.C/std::pow(iter,tp.t);

  while(!converged){

    for(i=0; i<numPar; i++){
      h.at(i) = (njm::rber(0.5) == 1 ? 1.0 : -1.0) * cm;
      parPH.at(i) = par.at(i) + h.at(i);
      parMH.at(i) = par.at(i) - h.at(i);
    }


    
    agent.tp.putPar(parPH);
    qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,agent);
    qEval.buildD1();
    try{
      qEval.solve();
    }
    catch(int e){
      std::cout << "Decomposition failed on thread "
		<< omp_get_thread_num() << "."
		<< std::endl
		<< "Lambda value is " << qEval.tp.lambda
		<< std::endl;
      throw(1);
    }
    valP = qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst,agent);

    agent.tp.putPar(parMH);
    qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,agent);
    qEval.buildD1();
    try{
      qEval.solve();
    }
    catch(int e){
      std::cout << "Decomposition failed on thread "
		<< omp_get_thread_num() << "."
		<< std::endl
		<< "Lambda value is " << qEval.tp.lambda
		<< std::endl;

      // save as much data as possible
      
      njm::toFile(njm::toString(qEval.D0,"\n",64,32),
		  njm::sett.datExt("D0.txt"));
      njm::toFile(njm::toString(qEval.D1,"\n",64,32),
		  njm::sett.datExt("D1.txt"));
      njm::toFile(njm::toString(qEval.D,"\n",64,32),
		  njm::sett.datExt("D.txt"));
      njm::toFile(njm::toString(qEval.R,"\n",64,32),
		  njm::sett.datExt("R.txt"));

      std::string name;
      for(int k = 0; k < qEval.numNodes; ++k){
	name = "D0L_" + njm::toString(k,"") + ".txt";
	njm::toFile(njm::toString(qEval.D0L.at(k),"\n",64,32),
		    njm::sett.datExt(name));

	name = "D1L_" + njm::toString(k,"") + ".txt";
	njm::toFile(njm::toString(qEval.D1L.at(k),"\n",64,32),
		    njm::sett.datExt(name));
	
	name = "RL_" + njm::toString(k,"") + ".txt";
	njm::toFile(njm::toString(qEval.RL.at(k),"\n",64,32),
		    njm::sett.datExt(name));
	
      	name = "phiL_" + njm::toString(k,"") + ".txt";
	njm::toFile(njm::toString(qEval.phiL.at(k),"\n",64,32),
		    njm::sett.datExt(name));
      }
      for(int t = 0; t < s.sD.time; ++t){
	for(int k = 0; k < qEval.numNodes; ++k){
	  name= "phiPsiTL_" + njm::toString(t,"") +
	    "_" + njm::toString(k,"") + ".txt";
	  njm::toFile(njm::toString(qEval.phiPsiTL.at(t).at(k),"\n",64,32),
		      njm::sett.datExt(name));
	}
      }
      
      throw(1);
    }
    valM = qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst,agent);

    
    for(i=0; i<numPar; i++)
      par.at(i) = par.at(i) - mu*(valP - valM)/(2.0*h.at(i));

    
    // if(omp_get_thread_num() == 0)
    //   std::cout << "iter: " + njm::toString(iter,"",4,0) +
    // 	" || " + njm::toString(valP,"",24,16) + " - " +
    // 	njm::toString(valM,"",24,16) + " -> " +
    // 	njm::toString(mu,"",6,4) + " , " + njm::toString(cm,"",6,4) +
    // 	" || " + njm::toString(par,", ","") << "\r" << std::flush;


    ++iter;
      
    mu = tp.A/std::pow(tp.B + iter,tp.ell);
    cm = tp.C/std::pow(iter,tp.t);

    
    if(mu < tp.muMin){
      converged = 1;
    }
    
  }

  agent.tp.putPar(par); // assign optimized par to the agent
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
  dim = (tp.dfLat * tp.dfLong + 1) * lenPsi;
  
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


  std::vector<double> mdsLat = fD.centroidsMdsLat;
  std::vector<double> mdsLong = fD.centroidsMdsLong;

  std::sort(mdsLat.begin(),mdsLat.end());
  std::sort(mdsLong.begin(),mdsLong.end());
  
  gsl_vector *knotsLat, *knotsLong;
  knotsLat = gsl_vector_alloc(nbreakLat);
  knotsLong = gsl_vector_alloc(nbreakLong);

  int ind;
  for(i = 0; i < nbreakLat; ++i){
    ind = (i*(fD.numNodes-1))/(nbreakLat-1);
    gsl_vector_set(knotsLat,i,mdsLat.at(ind));
  }
  for(i = 0; i < nbreakLong; ++i){
    ind = (i*(fD.numNodes-1))/(nbreakLong-1);
    gsl_vector_set(knotsLong,i,mdsLong.at(ind));
  }
  
  gsl_bspline_knots(knotsLat,bsLat);
  gsl_bspline_knots(knotsLong,bsLong);

  gsl_vector *bLat, *bLong;
  bLat = gsl_vector_alloc(tp.dfLat);
  bLong = gsl_vector_alloc(tp.dfLong);


  phiL.clear();
  int u,v;
  double bb;
  for(i = 0; i < fD.numNodes; ++i){
    Eigen::SparseMatrix<double> phi(dim,lenPsi);
    
    gsl_bspline_eval(fD.centroidsMdsLat.at(i),
		     bLat,bsLat);
    gsl_bspline_eval(fD.centroidsMdsLong.at(i),
		     bLong,bsLong);
    for(j = 0; j < lenPsi; ++j){
      phi.insert(j*(tp.dfLat*tp.dfLong + 1),j) = 1.0;
      for(u = 0; u < tp.dfLat; ++u){
	for(v = 0; v < tp.dfLong; ++v){
	  bb = gsl_vector_get(bLat,u)*gsl_vector_get(bLong,v);
	  if(bb != 0.0)
	    phi.insert(j*(tp.dfLat*tp.dfLong+1) + u*tp.dfLong + v + 1,j) = bb;
	}
      }
    }
    phiL.push_back(phi);
  }

  gsl_bspline_free(bsLat);
  gsl_bspline_free(bsLong);
  gsl_vector_free(bLat);
  gsl_vector_free(bLong);
  gsl_vector_free(knotsLat);
  gsl_vector_free(knotsLong);
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


    // using the current values build D0
    f.preCompData(sDt,tDt,fD,dD,m,mP);
    f.getFeatures(sDt,tDt,fD,dD,m,mP);
    features=feat2Vec(fD.numNodes,sDt.status);


    phiPsiL=featToPhiPsi(features,fD.numNodes);
    
    phiPsiTL.push_back(phiPsiL);

    
    for(i = 0; i < fD.numNodes; ++i){
      D0L.at(i) += phiPsiL.at(i) * phiPsiL.at(i).transpose();
    }

    numNewInfec = sDt.newInfec.size();
    for(i = 0; i < numNewInfec; ++i)
      RL.at(i) += phiPsiL.at(sDt.newInfec.at(i)) * (1.0/double(fD.numNodes));
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
      for(j=0; j<tp.polReps; j++){
	std::fill(tDt.a.begin(),tDt.a.end(),0);
	std::fill(tDt.p.begin(),tDt.p.end(),0);  
	a.applyTrt(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);

	f.preCompData(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	f.getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
	features = feat2Vec(fD.numNodes,sD1T.at(t).status);

	phiPsiL = featToPhiPsi(features,fD.numNodes);
	if(j == 0)
	  phiPsiLavg = phiPsiL;
	else
	  for(k = 0; k < fD.numNodes; ++k)
	    phiPsiLavg.at(k) += phiPsiL.at(k);
      }
      for(k = 0; k < fD.numNodes; ++k)
	phiPsiLavg.at(k) /= double(tp.polReps);
    }
    else{
      f.preCompData(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      f.getFeatures(sD1T.at(t),tDt,fD,dD1T.at(t),m,mP);
      features = feat2Vec(fD.numNodes,sD1T.at(t).status);

      phiPsiLavg = featToPhiPsi(features,fD.numNodes);
    }

    for(k = 0; k < fD.numNodes; ++k)
      D1L.at(k) += phiPsiTL.at(t).at(k) * phiPsiLavg.at(k).transpose();
  }


  for(k = 0; k < fD.numNodes; ++k)
    D1L.at(k) *= tp.gamma; // discount factor

}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
solve(){
  Eigen::SparseMatrix<double> P(dim,dim);

  // P.setIdentity();

  int i,ind = tp.dfLat*tp.dfLong + 1;
  for(i = 0; i < dim; ++i)
    if(i % ind != 0)
      P.insert(i,i) = 1.0;

  // Eigen::SparseMatrix<double> DD = D.transpose() * D;
  // Eigen::SparseQR<Eigen::SparseMatrix<double>,
  // 		  Eigen::COLAMDOrdering<int> > DDlu(DD);
  // std::cout << "rank: " << DDlu.rank() << std::endl;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  solver.compute(D.transpose() * D + tp.lambda*P);

  if(solver.info() != Eigen::Success){
    // std::cout << "In M2QEval::solve(): decomposition failed."
    // 	      << std::endl;
    throw(1);
  }

  beta = solver.solve(-D.transpose() * R);

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
  
  int i,I = nodes.size();
  R.resize(dim);
  D0.resize(dim,dim);
  D1.resize(dim,dim);

  R.setZero();
  D0.setZero();
  D1.setZero();
  for(i = 0; i < I; ++i)
    D0 += D0L.at(nodes.at(i));
  for(i = 0; i < I; ++i)
    D1 += D1L.at(nodes.at(i));
  for(i = 0; i < I; ++i)
    R += RL.at(nodes.at(i));

  D = D1 - D0;

}





template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
buildD1(){
  std::vector<int> nodes;
  nodes.resize(numNodes);
  int i = 0;
  std::for_each(nodes.begin(),nodes.end(),
		[&i](int & x){
		  x = i++;
		});
  
  buildD1(nodes);
}


template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
buildD1(const std::vector<int> nodes){
  
  int i,I = nodes.size();
  D1.resize(dim,dim);

  D1.setZero();
  for(i = 0; i < I; ++i)
    D1 += D1L.at(nodes.at(i));

  D = D1 - D0;

}



template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
getRD(Eigen::VectorXd & R,
      Eigen::SparseMatrix<double> & D){
  R = this->R;
  D = this->D;
}


template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
setRD(const Eigen::VectorXd & R,
      const Eigen::SparseMatrix<double> & D){
  this->R = R;
  this->D = D;
}


template <class S, class A, class F,
	  class M,class MP>
void M2QEval<S,A,F,M,MP>::
tune(const std::vector<int> & status){
  int i,b,bS = (tp.bootSize * numNodes + 1);


  // all the nodes in a vector
  std::vector<int> nodes;
  nodes.reserve(numNodes);
  for(i = 0; i < numNodes; ++i)
    nodes.push_back(i);

  // get the powers on lambda
  int numLambdaPows = 15;
  std::vector<double> lambdaPows;
  for(i = 0; i < numLambdaPows; ++i)
    lambdaPows.push_back(i+10);

  // setup CV error container
  std::vector<double> lambdaCV;
  lambdaCV.resize(numLambdaPows);
  std::fill(lambdaCV.begin(),lambdaCV.end(),0.0);

  // setup CV node containers
  std::vector<std::vector<int> > selNodesTrain;
  std::vector<std::vector<int> > selNodesTest;
  selNodesTrain.resize(tp.bootReps);
  selNodesTest.resize(tp.bootReps);


  // sample the CV nodes
  std::pair<double,int> top;
  int trainInf,testInf;
  for(b = 0; b < tp.bootReps; ++b){
    std::priority_queue<std::pair<double,int> > ordNodes;
    for(i = 0; i < numNodes; ++i)
      ordNodes.push(std::pair<double,int>(njm::runif01(),nodes.at(i)));

    trainInf = 0;
    testInf = 0;

    selNodesTrain.at(b).clear();
    selNodesTest.at(b).clear();
    for(i = 0; i < numNodes; ++i){
      top = ordNodes.top();
      ordNodes.pop();
      if(i < bS)
	selNodesTrain.at(b).push_back(top.second);
      else
	selNodesTest.at(b).push_back(top.second);

      if(status.at(top.second) >= 2){
	if(i < bS)
	  ++trainInf;
	else
	  ++testInf;
      }
    }
    
    // std::cout << "Inf: [" << double(trainInf)/double(bS)
    // 	      << ", " << double(testInf)/double(numNodes-bS)
    // 	      << "]" << std::endl;
  }
    


  // D,R containers for testing and training
  Eigen::VectorXd Rtrain,Rtest;
  Eigen::SparseMatrix<double> Dtrain,Dtest;

  // CV bellman error
  for(b = 0; b < tp.bootReps; ++b){

    buildRD(selNodesTrain.at(b));
    getRD(Rtrain,Dtrain);

    buildRD(selNodesTest.at(b));
    getRD(Rtest,Dtest);

    for(i = 0; i < numLambdaPows; ++i){
      tp.lambda = std::pow(2.0,lambdaPows.at(i));

      try{
	setRD(Rtrain,Dtrain);
	solve();
	
	setRD(Rtest,Dtest);
	lambdaCV.at(i) += bellRes();
      }
      catch(int e){
	lambdaCV.at(i) += std::numeric_limits<double>::max();
      }
      
      // std::cout << "    lambda (" << tp.lambda << ") -> "
      // 		<< lambdaCV.at(i) << std::endl;
    }
  }


  // pick off the best
  std::priority_queue<std::pair<double,double> > ordLambda;
  for(i = 0; i < numLambdaPows; ++i)
    ordLambda.push(std::pair<double,double>(-lambdaCV.at(i),lambdaPows.at(i)));
  double bestPow = ordLambda.top().second;

  tp.lambda = std::pow(2.0,bestPow);
  // std::cout << "lambda mid (" << tp.lambda << ") -> "
  // 	    << -ordLambda.top().first << std::endl;

  // now grid it up finer
  lambdaPows.clear();
  lambdaCV.clear();
  
  numLambdaPows = 10 + 1;
  for(i = 0; i < numLambdaPows; ++i)
    if(i != 5)
      lambdaPows.push_back(bestPow + double(i - 5)*0.4);
  --numLambdaPows;


  // error for new lambdas
  lambdaCV.resize(numLambdaPows);
  std::fill(lambdaCV.begin(),lambdaCV.end(),0.0);
  for(b = 0; b < tp.bootReps; ++b){

    buildRD(selNodesTrain.at(b));
    getRD(Rtrain,Dtrain);

    buildRD(selNodesTest.at(b));
    getRD(Rtest,Dtest);

    for(i = 0; i < numLambdaPows; ++i){
      tp.lambda = std::pow(2.0,lambdaPows.at(i));

      try{
	setRD(Rtrain,Dtrain);
	solve();
	
	setRD(Rtest,Dtest);
	lambdaCV.at(i) += bellRes();
      }
      catch(int e){
	lambdaCV.at(i) += std::numeric_limits<double>::max();
      }
      
      // std::cout << "    lambda (" << tp.lambda << ") -> "
      // 		<< lambdaCV.at(i) << std::endl;
    }
  }


  for(i = 0; i < numLambdaPows; ++i)
    ordLambda.push(std::pair<double,double>(-lambdaCV.at(i),lambdaPows.at(i)));
  tp.lambda = std::pow(2.0,ordLambda.top().second);
  
  // std::cout << "lambda final (" << tp.lambda << ") -> "
  // 	    << -ordLambda.top().first << std::endl;
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

  
  return (phiPsi.transpose()*beta).sum();
}
