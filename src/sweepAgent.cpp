#include "sweepAgent.hpp"


template class SweepAgent<GravityModel,GravityParam>;



template <class Model, class ModelParam>
SweepAgent<Model,ModelParam>::SweepAgent(){
  tp.sigma = 906.3819; // standard deviation of wnsD.txt
  tp.weights.ones(4);
  tp.maxSweep=10;
  name="sweep";
}



template <class Model, class ModelParam>
int SweepAgent<Model,ModelParam>::numFeatures = 4;



template<class Model, class ModelParam>
void SweepAgent<Model, ModelParam>::applyTrt(const SimData & sD,
					     TrtData & tD,
					     const FixedData & fD,
					     const DynamicData & dD,
					     const Model & m,
					     ModelParam & mP){
  if(sD.notInfec.empty())
    return;
  
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);
  
  m.load(sD,tD,fD,dD,mP);
  getFeatures(sD,tD,fD,dD,m,mP);

  double bestValue = calcValue(),curValue;

  int bestInd=0,i,j,first;

  std::vector<double> curNot,curInf;

  /////// obtain initial treatments
  // preventative
  for(i=0; i<numPre; i++){
    first=1;
    for(j=0; j<sD.numNotInfec; j++){
      if(first && !tD.p.at(sD.notInfec.at(j))){
	first=0;
	bestInd=j;
	updateSameP(j,1,sD,tD,fD,dD,m,mP);
	bestValue = calcValue();
	updateSameP(j,-1,sD,tD,fD,dD,m,mP);
      }
      else if(!tD.p.at(sD.notInfec.at(j))){
	updateSameP(j,1,sD,tD,fD,dD,m,mP);
	curValue = calcValue();
	updateSameP(j,-1,sD,tD,fD,dD,m,mP);
	if(curValue < bestValue){
	  bestInd=j;
	  bestValue=curValue;
	}
      }
    }

    curNot.push_back(bestInd);
    updateSameP(bestInd,1,sD,tD,fD,dD,m,mP);
  }

  for(i=0; i<numAct; i++){
    first=1;
    for(j=0; j<sD.numInfected; j++){
      if(first && !tD.a.at(sD.infected.at(j))){
	first=0;
	bestInd=j;
	updateSameA(j,1,sD,tD,fD,dD,m,mP);
	bestValue = calcValue();
	updateSameA(j,-1,sD,tD,fD,dD,m,mP);
      }
      else if(!tD.a.at(sD.infected.at(j))){
	updateSameA(j,1,sD,tD,fD,dD,m,mP);
	curValue = calcValue();
	updateSameA(j,-1,sD,tD,fD,dD,m,mP);
	if(curValue < bestValue){
	  bestInd=j;
	  bestValue=curValue;
	}
      }
    }

    curInf.push_back(bestInd);
    updateSameA(bestInd,1,sD,tD,fD,dD,m,mP);
  }


  
  
  //////// start the sweep
  int curSweep=0,change=1,lastInd;
  while(change && curSweep++ < tp.maxSweep){
    change=0;
    
    // sweep preventative
    for(i=0; i<numPre; i++){
      bestInd = lastInd = curNot.at(i);
      for(j=0; j<sD.numNotInfec; j++){
	if(!tD.p.at(sD.notInfec.at(j))){
	  updateSwapP(lastInd,j,sD,tD,fD,dD,m,mP);
	  lastInd = j;
	  curValue = calcValue();
	  if(curValue < bestValue){
	    bestInd = j;
	    bestValue = curValue;
	    change=1;
	  }
	}
      }

      curNot.at(i) = bestInd;
      updateSwapP(lastInd,bestInd,sD,tD,fD,dD,m,mP);
    }


    
    // sweep active
    for(i=0; i<numAct; i++){
      bestInd = lastInd = curInf.at(i);
      for(j=0; j<sD.numInfected; j++){
	if(!tD.a.at(sD.infected.at(j))){
	  updateSwapA(lastInd,j,sD,tD,fD,dD,m,mP);
	  lastInd = j;
	  curValue = calcValue();
	  if(curValue < bestValue){
	    bestInd = j;
	    bestValue = curValue;
	    change=1;
	  }
	}
      }

      curInf.at(i) = bestInd;
      updateSwapA(lastInd,bestInd,sD,tD,fD,dD,m,mP);
    }
  }
}



template<class Model, class ModelParam>
void SweepAgent<Model, ModelParam>::getFeatures(const SimData & sD,
						const TrtData & tD,
						const FixedData & fD,
						const DynamicData & dD,
						const Model & m,
						const ModelParam & mP){
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);

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
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  // double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD0)*fD.numNodes + *itD1)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // notFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      if(!tD.a.at(*itD1))
	p.push(-fD.dist.at((*itD0)*fD.numNodes + *itD1));
    if(p.empty())
      notFeat(i,featNum) = 0;
    else
      notFeat(i,featNum) = -p.top()/tp.sigma;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD1)*fD.numNodes + *itD0)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // infFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      if(!tD.p.at(*itD1))
	p.push(-fD.dist.at((*itD1)*fD.numNodes + *itD0));
    if(p.empty())
      infFeat(i,featNum) = 0;
    else
    infFeat(i,featNum) = -p.top()/tp.sigma;
  }

  featNum++;
}
							  


template<class Model, class ModelParam>
double SweepAgent<Model, ModelParam>::calcValue(){
  arma::colvec infValue,notValue;
  infValue = infFeat*tp.weights;
  notValue = notFeat*tp.weights;
  return arma::sum(infValue) + arma::sum(notValue);
}



template<class Model, class ModelParam>
void SweepAgent<Model, ModelParam>::updateSameP(const int ind, const int change,
						const SimData & sD,
						TrtData & tD,
						const FixedData & fD,
						const DynamicData & dD,
						const Model & m,
						ModelParam & mP){
  // update existing data
  if(change == 1)
    tD.p.at(sD.notInfec.at(ind)) = 1;
  else
    tD.p.at(sD.notInfec.at(ind)) = 0;
  mP.infProbsBase.col(ind) -= change*mP.trtPre;
  mP.infProbsSep.col(ind) = 1/(1+arma::exp(mP.infProbsBase.col(ind)));

  //////// update features
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;  
  // feature 0
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat(ind,featNum) = 1 - arma::prod(mP.infProbsSep.col(ind));

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
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  // double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD0)*fD.numNodes + *itD1)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // notFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      if(!tD.a.at(*itD1))
	p.push(-fD.dist.at((*itD0)*fD.numNodes + *itD1));
    if(p.empty())
      notFeat(i,featNum) = 0;
    else
      notFeat(i,featNum) = -p.top()/tp.sigma;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD1)*fD.numNodes + *itD0)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // infFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      if(!tD.p.at(*itD1))
	p.push(-fD.dist.at((*itD1)*fD.numNodes + *itD0));
    if(p.empty())
      infFeat(i,featNum) = 0;
    else
      infFeat(i,featNum) = -p.top()/tp.sigma;
  }

  featNum++;
}
			    


template<class Model, class ModelParam>
void SweepAgent<Model, ModelParam>::updateSwapP(const int oldInd,
						const int newInd,
						const SimData & sD,
						TrtData & tD,
						const FixedData & fD,
						const DynamicData & dD,
						const Model & m,
						ModelParam & mP){
  // update existing data
  tD.p.at(sD.notInfec.at(oldInd)) = 0;
  mP.infProbsBase.col(oldInd) += mP.trtPre;
  mP.infProbsSep.col(oldInd) = 1/(1+arma::exp(mP.infProbsBase.col(oldInd)));

  tD.p.at(sD.notInfec.at(newInd)) = 1;  
  mP.infProbsBase.col(newInd) -= mP.trtPre;
  mP.infProbsSep.col(newInd) = 1/(1+arma::exp(mP.infProbsBase.col(newInd)));

  // update features
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;    
  // feature 0
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat(oldInd,featNum) = 1 - arma::prod(mP.infProbsSep.col(oldInd));
  notFeat(oldInd,featNum) = 1 - arma::prod(mP.infProbsSep.col(newInd));  

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
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  // double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD0)*fD.numNodes + *itD1)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // notFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      if(!tD.a.at(*itD1))
	p.push(-fD.dist.at((*itD0)*fD.numNodes + *itD1));
    if(p.empty())
      notFeat(i,featNum) = 0;
    else
      notFeat(i,featNum) = -p.top()/tp.sigma;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD1)*fD.numNodes + *itD0)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // infFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      if(!tD.p.at(*itD1))
	p.push(-fD.dist.at((*itD1)*fD.numNodes + *itD0));
    if(p.empty())
      infFeat(i,featNum) = 0;
    else
      infFeat(i,featNum) = -p.top()/tp.sigma;
  }

  featNum++;
}




template<class Model, class ModelParam>
void SweepAgent<Model, ModelParam>::updateSameA(const int ind, const int change,
						const SimData & sD,
						TrtData & tD,
						const FixedData & fD,
						const DynamicData & dD,
						const Model & m,
						ModelParam & mP){
  // update existing data
  if(change == 1)
    tD.a.at(sD.infected.at(ind)) = 1;
  else
    tD.a.at(sD.infected.at(ind)) = 0;
  mP.infProbsBase.row(ind) -= change*mP.trtAct;
  mP.infProbsSep.row(ind) = 1/(1+arma::exp(mP.infProbsBase.row(ind)));

  // update features
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;    
  // feature 0
  infFeat(ind,featNum) = 1 - arma::prod(mP.infProbsSep.row(ind));  
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
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  // double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD0)*fD.numNodes + *itD1)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // notFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      if(!tD.a.at(*itD1))
	p.push(-fD.dist.at((*itD0)*fD.numNodes + *itD1));
    if(p.empty())
      notFeat(i,featNum) = 0;
    else
      notFeat(i,featNum) = -p.top()/tp.sigma;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD1)*fD.numNodes + *itD0)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // infFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      if(!tD.p.at(*itD1))
	p.push(-fD.dist.at((*itD1)*fD.numNodes + *itD0));
    if(p.empty())
      infFeat(i,featNum) = 0;
    else
      infFeat(i,featNum) = -p.top()/tp.sigma;
  }

  featNum++;
}
			    


template<class Model, class ModelParam>
void SweepAgent<Model, ModelParam>::updateSwapA(const int oldInd,
						const int newInd,
						const SimData & sD,
						TrtData & tD,
						const FixedData & fD,
						const DynamicData & dD,
						const Model & m,
						ModelParam & mP){
  // update existing data
  tD.a.at(sD.infected.at(oldInd)) = 0;
  mP.infProbsBase.row(oldInd) += mP.trtAct;
  mP.infProbsSep.row(oldInd) = 1/(1+arma::exp(mP.infProbsBase.row(oldInd)));

  tD.a.at(sD.infected.at(newInd)) = 1;  
  mP.infProbsBase.row(newInd) -= mP.trtAct;
  mP.infProbsSep.row(newInd) = 1/(1+arma::exp(mP.infProbsBase.row(newInd)));

  // update features
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;    
  // feature 0
  infFeat(oldInd,featNum) = 1 - arma::prod(mP.infProbsSep.row(oldInd));
  infFeat(newInd,featNum) = 1 - arma::prod(mP.infProbsSep.row(newInd));  
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
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  // double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD0)*fD.numNodes + *itD1)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // notFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      if(!tD.a.at(*itD1))
	p.push(-fD.dist.at((*itD0)*fD.numNodes + *itD1));
    if(p.empty())
      notFeat(i,featNum) = 0;
    else
      notFeat(i,featNum) = -p.top()/tp.sigma;
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    // totalDist=0;
    // for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
    //   totalDist += std::exp(fD.dist.at((*itD1)*fD.numNodes + *itD0)
    // 			    /tp.sigma);
    // totalDist /= fD.numNodes*tp.sigma;
    // infFeat(i,featNum) = totalDist;
    p=std::priority_queue<double>();
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      if(!tD.p.at(*itD1))
	p.push(-fD.dist.at((*itD1)*fD.numNodes + *itD0));
    if(p.empty())
      infFeat(i,featNum) = 0;
    else
      infFeat(i,featNum) = -p.top()/tp.sigma;
  }

  featNum++;
}











std::vector<double> SweepTuneParam::getPar() const {
  std::vector<double> par;
  par = arma::conv_to< std::vector<double> >::from(weights);
  // par.push_back(sigma);
  return par;
}


void SweepTuneParam::putPar(const std::vector<double> & par){
  // sigma = par.back();
  weights = arma::conv_to<arma::colvec>::from(par);
  // weights.resize(weights.n_elem - 1);
}
