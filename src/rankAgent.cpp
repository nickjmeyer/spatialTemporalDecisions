#include "rankAgent.hpp"


template class RankAgent<ToyFeatures0<GravityModel,GravityParam>,
			 GravityModel,GravityParam>;


template < class F, class M, class MP>
RankAgent<F,M,MP>::RankAgent(){
  tp.weights.ones(4);
  tp.numChunks = 3;
  tp.valReps=10;
  name="rank";
}

  
template < class F, class M, class MP>
int RankAgent<F,M,MP>::numFeatures = 4;


template < class F, class M, class MP>
void RankAgent<F,M,MP>::applyTrt(const SimData & sD,
				 TrtData & tD,
				 const FixedData & fD,
				 const DynamicData & dD,
				 const M & m,
				 MP & mP){
  if(sD.notInfec.empty())
    return;
  
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);
  
  m.load(sD,tD,fD,dD,mP);
  getFeatures(sD,tD,fD,dD,m,mP);

  std::priority_queue<std::pair<double,int> > sortInfected,sortNotInfec;
  
  int i,j,node0,addPre,addAct;
  int cI=0,cN=0;
  for(i=0; i<tp.numChunks; i++){

    infRanks = infFeat * tp.weights;
    notRanks = notFeat * tp.weights;


    for(j=0; j<sD.numInfected; j++){
      if(tD.a.at(sD.infected.at(j)))
	sortInfected.push(std::pair<double,int>(std::numeric_limits<double>
						::min(),j));
      else
	sortInfected.push(std::pair<double,int>(infRanks(j),j));
    }
    for(j=0; j<sD.numNotInfec; j++){
      if(tD.p.at(sD.notInfec.at(j)))
	sortNotInfec.push(std::pair<double,int>(std::numeric_limits<double>
						::min(),j));
      else
	sortNotInfec.push(std::pair<double,int>(notRanks(j),j));
    }


    addAct = (int)((i+1)*numAct/std::min(tp.numChunks,numAct)) -
      (int)(i*numAct/std::min(tp.numChunks,numAct));
    for(; cI<(cI+addAct) && cI<numAct; cI++){
      node0=sortInfected.top().second;
      tD.a.at(sD.infected.at(node0)) = 1;
      mP.infProbsBase.row(node0) -= mP.trtAct;
      mP.infProbsSep.row(node0) = 1/(1+arma::exp(mP.infProbsBase.row(node0)));
      sortInfected.pop();
      
    }
    
    addPre = (int)((i+1)*numPre/std::min(tp.numChunks,numPre)) -
	     (int)(i*numPre/std::min(tp.numChunks,numPre)); 
    for(; cN<(cN+addPre) && cN<numPre; cN++){
      node0=sortNotInfec.top().second;
      tD.p.at(sD.notInfec.at(node0)) = 1;
      mP.infProbsBase.col(node0) -= mP.trtPre;
      mP.infProbsSep.col(node0) = 1/(1+arma::exp(mP.infProbsBase.col(node0)));
      sortNotInfec.pop();
    }

    
    if((i+1) < tp.numChunks){
      
      getFeatures(sD,tD,fD,dD,m,mP);
    }
  }
}




template < class F, class M, class MP>
void RankAgent<F,M,MP>::getFeatures(const SimData & sD,
				    const TrtData & tD,
				    const FixedData & fD,
				    const DynamicData & dD,
				    const M & m,
				    const MP & mP){
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);

  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1;


  // feature 0
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,0).t();
  
  featNum++;
  
  // feature 1
  SystemLight<M,MP> s(sD,tD,fD,dD,m,mP);
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

  double minDist,distDepth;
  for(i=0,itD0=sD.notInfec.begin(); i<sD.numNotInfec; i++,itD0++){
    minDist=0;
    for(j=0,itD1=sD.infected.begin(); j<sD.numInfected; j++,itD1++)
      if((1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1))) > minDist)
	minDist = 1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1));
    distDepth = minDist/(halfPlaneDepth(fD.centroidsLong.at(*itD0),
					fD.centroidsLat.at(*itD0),
					numInfNoTrt,
					infNoTrtLong,
					infNoTrtLat)
			 + (1.0/(double)(numInfNoTrt+1)));
    notFeat(i,featNum) = std::log(1.0+distDepth);
  }
  for(i=0,itD0=sD.infected.begin(); i<sD.numInfected; i++,itD0++){
    minDist=0;
    for(j=0,itD1=sD.notInfec.begin(); j<sD.numNotInfec; j++,itD1++)
      if((1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1))) > minDist)
	minDist = 1.0/fD.logDist.at((*itD0)*fD.numNodes+(*itD1));
    distDepth = minDist/(halfPlaneDepth(fD.centroidsLong.at(*itD0),
					fD.centroidsLat.at(*itD0),
					numNotNoTrt,
					notNoTrtLong,
					notNoTrtLat)
			 + (1.0/(double)(numNotNoTrt+1)));
    infFeat(i,featNum) = std::log(1.0+distDepth);
  }
  featNum++;
  

  // feature 3
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
    
}






std::vector<double> RankTuneParam::getPar() const {
  std::vector<double> par;
  par = arma::conv_to< std::vector<double> >::from(weights);
  // par.push_back(sigma);
  return par;
}


void RankTuneParam::putPar(const std::vector<double> & par){
  // sigma = par.back();
  weights = arma::conv_to<arma::colvec>::from(par);
  // weights.resize(weights.n_elem - 1);
}
