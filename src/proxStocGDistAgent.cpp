#include "proxStocGDistAgent.hpp"

template class ProxStocGDistAgent<ModelGravityGDist>;

template class ProxStocGDistAgent<ModelTimeGDistTrendPow>;

template class ProxStocGDistAgent<ModelTimeExpCavesGDistTrendPowCon>;

template class ProxStocGDistAgent<ModelTimeExpCavesGDist>;

template class ProxStocGDistAgent<ModelRadius>;

template class ProxStocGDistAgent<ModelGDistKern>;


std::vector<double> ProxStocGDistTuneParam::getPar() const {
  return std::vector<double>();
}


void ProxStocGDistTuneParam::putPar(const std::vector<double> & par){
}


template<class M>
ProxStocGDistAgent<M>::ProxStocGDistAgent(){
  tp.corrGoal = 0.8;
}


template<class M>
std::string ProxStocGDistAgent<M>::name = "proxStoc";

template<class M>
void ProxStocGDistAgent<M>::reset(){
}

template <class M>
void ProxStocGDistAgent<M>::applyTrt(const SimData & sD,
				TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD,
				M & m){
  numPre = getNumPre(sD,tD,fD,dD);
  numAct = getNumAct(sD,tD,fD,dD);

  int i,j,node0,node1;
  double minDist,curDist,maxDist;

  maxDist = std::numeric_limits<double>::max();

  std::vector<double> notInfecVals,infectedVals;
  for(i=0; i<sD.numNotInfec; i++){
    node0 = sD.notInfec.at(i);

    minDist=maxDist;
    for(j=0; j<sD.numInfected; j++){
      node1=sD.infected.at(j);
      curDist=fD.gDist.at(node0*fD.numNodes + node1);
      if(minDist > curDist)
	minDist = curDist;
    }

    notInfecVals.push_back(minDist);
  }

  
  for(i=0; i<sD.numInfected; i++){
    node0=sD.infected.at(i);

    minDist=maxDist;
    for(j=0; j<sD.numNotInfec; j++){
      node1=sD.notInfec.at(j);
      curDist=fD.gDist.at(node0*fD.numNodes + node1);
      if(minDist > curDist)
	minDist = curDist;
    }

    infectedVals.push_back(minDist);
  }

  // variance calculation
  double varNot = 1.0, varInf = 1.0;
  if(sD.numNotInfec > 1)
    varNot = optCorr(notInfecVals,tp.corrGoal);

  if(sD.numInfected > 1)
    varInf = optCorr(infectedVals,tp.corrGoal);

  std::vector<double> notInfecScore = getScores(notInfecVals);
  std::vector<double> infectedScore = getScores(infectedVals);

  std::priority_queue<std::pair<double,int> > sortInfected;
  std::priority_queue<std::pair<double,int> > sortNotInfec;

  for(i = 0; i < sD.numInfected; ++i){
    node0 = sD.infected.at(i);
    double val = varInf*njm::rnorm01() + infectedScore.at(i);
    sortInfected.push(std::pair<double,int>(val,node0));
  }
	 
  for(i = 0; i < sD.numNotInfec; ++i){
    node0 = sD.notInfec.at(i);
    double val = varNot*njm::rnorm01() + notInfecScore.at(i);
    sortNotInfec.push(std::pair<double,int>(val,node0));
  }

  for(i=0; i<numAct; i++){
    tD.a.at(sortInfected.top().second) = 1;
    sortInfected.pop();
  }

  for(i=0; i<numPre; i++){
    tD.p.at(sortNotInfec.top().second) = 1;
    sortNotInfec.pop();
  }
}



double optCorr(const std::vector<double> vals,
	       const double goal){
  int i,I = 1000,j,J = int(vals.size());

  // generate std normals
  std::vector<std::vector<double> > stdNorm(I);
  for(i = 0; i < I; ++i){
    stdNorm.at(i).clear();
    for(j = 0; j < J; ++j){
      stdNorm.at(i).push_back(njm::rnorm01());
    }
  }

  
  double var = 1.0;
  double corr = getCorr(vals,stdNorm);

  double mod = 1.0;
  double dec = 0.99;
  bool pOver = corr > goal, over;
  unsigned int iter = 0, maxIter = 1000;
  while(std::abs(goal-corr) > 0.001 && iter++ < maxIter){
    over = corr > goal;
    if(over ^ pOver)
      mod *= dec;

    if(over)
      var *= 1.0 + mod;
    else
      var *= 1.0/(1.0 + mod);

    pOver = over;


    corr = getCorr(vals,appScale(stdNorm,var));
  }
  
  return var;
}


std::vector<double> getScores(std::vector<double> vals){
  std::for_each(vals.begin(),vals.end(),
		[](double & x){
		  return -std::log(x+1.0);
		});
  return vals;
}


std::vector<std::vector<double> >
appScale(std::vector<std::vector<double> > vals,
	 const double scale){
  std::for_each(vals.begin(),vals.end(),
		[scale](std::vector<double> & x){
		  std::for_each(x.begin(),x.end(),
				[scale](double & x){
				  x *= scale;
				});
		});
  return vals;
}


double getCorr(const std::vector<double> vals,
	       const std::vector<std::vector<double> > jitter){
  std::vector<double> scores = getScores(vals), scoresAdd;
  int i,I = int(jitter.size());
  int j,J = int(vals.size());

  double corr = 0.0;
  for(i = 0; i < I; ++i){
    scoresAdd = scores;
    for(j = 0; j < J; ++j){
      scoresAdd.at(j) += jitter.at(i).at(j);
    }

    corr += getCorr(getRank(scores),getRank(scoresAdd));
  }
  return corr/double(I);
}


double getCorr(const std::vector<double> vals0,
	       const std::vector<double> vals1){
  int i,I = int(vals0.size());
  RunStats rs0(vals0),rs1(vals1);
  double corr = 0.0;
  for(i = 0; i < I; ++i){
    double v0,v1;
    v0 = (vals0.at(i) - rs0.smean())/rs0.ssd();
    v1 = (vals1.at(i) - rs1.smean())/rs1.ssd();
    corr += v0*v1;
  }
  return corr / double(I-1);
}


std::vector<double> getRank(const std::vector<double> scores){
  std::priority_queue<std::pair<double,int> > pq;
  int i,I = int(scores.size());
  for(i = 0; i < I; ++i){
    pq.push(std::pair<double,int>(scores.at(i),i));
  }

  std::vector<double> rank(I);
  for(i = 0; i < I; ++i){
    rank.at(pq.top().second) = i;
  }
  return rank;
}
