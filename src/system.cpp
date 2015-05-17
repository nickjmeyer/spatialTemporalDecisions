#include "system.hpp"

template class System<ModelGravityGDist,
		      ModelGravityGDist>;

template class System<ModelGravityGDistTrend,
		      ModelGravityGDistTrend>;

template class System<ModelGravityGDistTrendPow,
		      ModelGravityGDistTrendPow>;

template class System<ModelGravityGDistTrendPowCon,
		      ModelGravityGDistTrendPowCon>;

template class System<ModelGravityEDist,
		      ModelGravityEDist>;

template class System<ModelTimeGDist,
		      ModelTimeGDist>;

template class System<ModelTimeGDistTrend,
		      ModelTimeGDistTrend>;

template class System<ModelTimeGDistTrendPow,
		      ModelTimeGDistTrendPow>;

template class System<ModelTimeGDistTrendPowCon,
		      ModelTimeGDistTrendPowCon>;

template class System<ModelTimeExpGDist,
		      ModelTimeExpGDist>;

template class System<ModelTimeExpGDistTrend,
		      ModelTimeExpGDistTrend>;

template class System<ModelTimeExpGDistTrendPow,
		      ModelTimeExpGDistTrendPow>;

template class System<ModelTimeExpGDistTrendPowCon,
		      ModelTimeExpGDistTrendPowCon>;

template class System<ModelTimeExpCavesGDist,
		      ModelTimeExpCavesGDist>;

template class System<ModelTimeExpCavesGDistTrend,
		      ModelTimeExpCavesGDistTrend>;

template class System<ModelTimeExpCavesGDistTrendPow,
		      ModelTimeExpCavesGDistTrendPow>;

template class System<ModelTimeExpCavesGDistTrendPowCon,
		      ModelTimeExpCavesGDistTrendPowCon>;

template class System<ModelTimeExpCavesGDist,
		      ModelTimeExpCavesEDist>;

template class System<ModelTimeExpCavesEDist,
		      ModelTimeExpCavesEDist>;

template class System<ModelTimeExpCavesGDist,
		      ModelTimeGDist>;

template class System<ModelTimeExpCavesGDist,
		      ModelGravityGDist>;

template class System<ModelTimeExpCavesGDist,
		      ModelRadius>;

template class System<ModelRadius,
		      ModelRadius>;

template class System<ModelTimeExpCavesGDist,
		      ModelGDist>;

template class System<ModelTimeExpCavesGDist,
		      ModelCovar>;

template class System<ModelTimeExpCavesGDist,
		      ModelGDistKern>;

template class System<ModelGDist,
		      ModelGDist>;

template class System<ModelGDistPow,
		      ModelGDistPow>;

template class System<ModelGDistKern,
		      ModelGDistKern>;

template class System<ModelCovar,
		      ModelCovar>;

template class System<ModelRad,
		      ModelRad>;



template <class MG,
	  class ME>
System<MG,
       ME>::System(){
  specialInit = 0;
  
  initialize();
}


template <class MG,
	  class ME>
System<MG,
       ME>::System(const SimData & sD,
		   const TrtData & tD,
		   const FixedData & fD,
		   const DynamicData & dD,
		   const MG & modelGen,
		   const ME & modelEst){
  specialInit = 1;
    
  this->sD_r = sD;
  this->tD_r = tD;
  this->fD = fD;
  this->dD_r = dD;
  this->modelGen_r = modelGen;
  this->modelEst_r = modelEst;
  revert();
}


template <class MG,
	  class ME>
System<MG,
       ME>::System(const std::string file){
  specialInit = 1;
  
  initialize();

  std::vector<int> historyFile;
  njm::fromFile(historyFile,njm::sett.srcExt(file));


  int size = int(historyFile.size());
  int numPoints = size/fD.numNodes;

  // history
  sD_r.history.clear();
  sD_r.history.resize(numPoints - 1);
  
  int i,j,k;
  for(i = 0, k = 0; i < (numPoints - 1); i++){
    sD_r.history.at(i).clear();
    for(j = 0; j < fD.numNodes; j++, k++){
      sD_r.history.at(i).push_back(historyFile.at(k));
    }
  }


  // status
  sD_r.status.clear();
  for(j = 0; j < fD.numNodes; j++, k++)
    sD_r.status.push_back(historyFile.at(k));


  // current treatments
  tD_r.a.resize(fD.numNodes);
  tD_r.p.resize(fD.numNodes);
  for(j = 0; j < fD.numNodes; j++){
    tD_r.a.at(j) = tD_r.p.at(j) = 0;
    if(sD_r.status.at(j) == 1)
      tD_r.p.at(j) = 1;
    else if(sD_r.status.at(j) == 3)
      tD_r.a.at(j) = 1;
  }

  // past treatments
  tD_r.aPast.resize(fD.numNodes);
  tD_r.pPast.resize(fD.numNodes);
  if(numPoints > 1)
    for(j = 0; j < fD.numNodes; j++){
      tD_r.aPast.at(j) = tD_r.pPast.at(j) = 0;
      if(sD_r.history.at(numPoints - 2).at(j) == 1)
	tD_r.pPast.at(j) = 1;
      else if(sD_r.history.at(numPoints - 2).at(j) == 3)
	tD_r.aPast.at(j) = 1;
    }
  else{
    std::fill(tD_r.pPast.begin(),tD_r.pPast.end(),0);
    std::fill(tD_r.aPast.begin(),tD_r.aPast.end(),0);
  }
	
  // infected & not infected
  sD_r.infected.clear();
  sD_r.notInfec.clear();
  for(j = 0; j < fD.numNodes; j++)
    if(sD_r.status.at(j) < 2)
      sD_r.notInfec.push_back(j);
    else
      sD_r.infected.push_back(j);
  sD_r.numInfected = sD_r.infected.size();
  sD_r.numNotInfec = sD_r.notInfec.size();

  // newly infected
  sD_r.newInfec.clear();
  if(numPoints > 1){
    for(j = 0; j < fD.numNodes; j++)
      if(sD_r.status.at(j) >= 2 && sD_r.history.at(numPoints - 2).at(j) < 2)
	sD_r.newInfec.push_back(j);
  }
  else{
    sD_r.newInfec = sD_r.infected;
  }


  // time infected
  sD_r.timeInf.resize(fD.numNodes);
  std::fill(sD_r.timeInf.begin(),sD_r.timeInf.end(),0);
  for(i = 0; i < (numPoints - 1); i++)
    for(j = 0; j < fD.numNodes; j++)
      if(sD_r.history.at(i).at(j) >= 2)
	sD_r.timeInf.at(j)++;
  for(j = 0; j < fD.numNodes; j++)
    if(sD_r.status.at(j) >= 2)
      sD_r.timeInf.at(j)++;
  
  // current time
  sD_r.time = numPoints - 1;

  // load probs
  modelGen_r.setFill(sD_r,tD_r,fD,dD_r);

  revert();
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::reset(const std::vector<int> & ind){
  if(!specialInit){
    // reset SimData
    sD_r.time = 0;
    sD_r.numInfected = ind.size();
    sD_r.numNotInfec = fD.numNodes - sD_r.numInfected;

    sD_r.infected = ind;
    int i,j;
    sD_r.notInfec.clear();
    for(i = 0,j = 0; i < fD.numNodes; ++i){
      if(j < sD_r.numInfected && i == sD_r.infected.at(j))
	++j;
      else
	sD_r.notInfec.push_back(i);
    }

    sD_r.newInfec = sD_r.infected;
    sD_r.timeInf.resize(fD.numNodes);
    std::fill(sD_r.timeInf.begin(),sD_r.timeInf.end(),0);
    for(i = 0; i < sD_r.numInfected; ++i)
      sD_r.timeInf.at(sD_r.infected.at(i)) = 1;

    sD_r.status.resize(fD.numNodes);
    std::fill(sD_r.status.begin(),sD_r.status.end(),0);
    for(i = 0; i < sD_r.numInfected; ++i)
      sD_r.status.at(sD_r.infected.at(i)) = 2;

    sD_r.history.clear();


    // reset TrtData
    tD_r.a.resize(fD.numNodes);
    tD_r.p.resize(fD.numNodes);
    tD_r.aPast.resize(fD.numNodes);
    tD_r.pPast.resize(fD.numNodes);

    std::fill(tD_r.a.begin(),tD_r.a.end(),0);
    std::fill(tD_r.p.begin(),tD_r.p.end(),0);
    std::fill(tD_r.aPast.begin(),tD_r.aPast.end(),0);
    std::fill(tD_r.pPast.begin(),tD_r.pPast.end(),0);


    // reset DynamicData
    // nothing to do for this....DynamicData isn't used

    // load probs
    modelGen_r.setFill(sD_r,tD_r,fD,dD_r);
    
  }
  
  // revert
  revert();
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::revert(){
  sD = sD_r;
  tD = tD_r;
  dD = dD_r;

  modelGen = modelGen_r;
  modelEst = modelEst_r;
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::checkPoint(){
  sD_r = sD;
  tD_r = tD;
  dD_r = dD;
  
  modelGen_r = modelGen;
  modelEst_r = modelEst;
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::initialize(){
  njm::fromFile(fD.fips,njm::sett.srcExt("fips.txt"));
  fD.numNodes = fD.fips.size();
  njm::fromFile(fD.gDist,njm::sett.srcExt("gDist.txt"));
  njm::fromFile(fD.caves,njm::sett.srcExt("caves.txt"));
  njm::fromFile(fD.covar,njm::sett.srcExt("xcov.txt"));
  fD.numCovar = ((int)fD.covar.size())/fD.numNodes;
  njm::fromFile(fD.network,njm::sett.srcExt("network.txt"));
  
  njm::fromFile(fD.centroidsLong,njm::sett.srcExt("centroidsLong.txt"));
  njm::fromFile(fD.centroidsLat,njm::sett.srcExt("centroidsLat.txt"));
  njm::fromFile(fD.centroidsMdsLong,njm::sett.srcExt("centroidsMdsLong.txt"));
  njm::fromFile(fD.centroidsMdsLat,njm::sett.srcExt("centroidsMdsLat.txt"));
  
  njm::fromFile(fD.subGraph,njm::sett.srcExt("subGraph.txt"));
  njm::fromFile(fD.betweenness,njm::sett.srcExt("betweenness.txt"));

  njm::fromFile(fD.priorTrtMean,njm::sett.srcExt("priorTrtMean.txt"));

  // only start treatment at time trtStart and on  
  njm::fromFile(fD.trtStart,njm::sett.srcExt("trtStart.txt"));
  // only update every period steps
  njm::fromFile(fD.period,njm::sett.srcExt("period.txt"));
  // final time step in simulation
  njm::fromFile(fD.finalT,njm::sett.srcExt("finalT.txt"));

  preCompData();

  
  modelGen_r = MG(fD);
  modelGen_r.read();
  
  modelEst_r = ME(fD);
  
  modelGen_r.setType(INVALID);
  modelEst_r.setType(INVALID);

}


template <class MG,
	  class ME>
void System<MG,
	    ME>::preCompData(){
  int i,j,k,tot;

  // eDist
  double iLat,iLong,jLat,jLong;
  fD.eDist.clear();
  fD.eDist.reserve(fD.numNodes);
  for(i = 0; i < fD.numNodes; ++i){
    iLat = fD.centroidsLat.at(i);
    iLong = fD.centroidsLong.at(i);
    for(j = 0; j < fD.numNodes; ++j){
      jLat = fD.centroidsLat.at(j);
      jLong = fD.centroidsLong.at(j);
      fD.eDist.push_back(std::sqrt((iLat - jLat)*(iLat - jLat)
				   + (iLong - jLong)*(iLong - jLong)));
    }
  }

  // circle mass
  fD.cm.clear();
  fD.cm.reserve(fD.numNodes*fD.numNodes);
  for(i = 0; i < fD.numNodes; ++i){
    for(j = 0; j < fD.numNodes; ++j){
      double numCaves = 0;
      double d = fD.gDist[i*fD.numNodes + j];
      for(k = 0; k < fD.numNodes; ++k){
	if(fD.gDist[i*fD.numNodes + k] <= d){
	  numCaves += fD.caves[k];
	}
      }
      fD.cm.push_back(numCaves);
    }
  }

  double maxVal = std::numeric_limits<double>::lowest();

  // proportion of caves, (caves[i] + 1)/(max(caves) + 1)
  for(i = 0 ; i < fD.numNodes; ++i){
    if(fD.caves.at(i) > maxVal)
      maxVal = fD.caves.at(i);
  }
  fD.propCaves = fD.caves;
  std::for_each(fD.propCaves.begin(),fD.propCaves.end(),
		[&maxVal](double & x){
		  x=(x+1.0)/(maxVal+1.0);
		});

  // proprtion of log caves, (log(caves[i]+1)+1)/(log(max(caves)+1)+1)
  fD.logPropCaves = fD.caves;
  std::for_each(fD.logPropCaves.begin(),fD.logPropCaves.end(),
		[&maxVal](double & x){
		  x=(std::log(x+1.0)+1.0)/(std::log(maxVal+1.0)+1.0);
		});

  // rank of caves
  int numGt;
  double curCaves;
  fD.rankCaves = std::vector<double> (fD.numNodes,0);
  for(i = 0; i < fD.numNodes; ++i){
    curCaves = fD.caves.at(i);
    numGt = 0;
    for(j = 0; j < fD.numNodes; ++j)
      if(fD.caves.at(j) >= curCaves)
	++numGt;
    fD.rankCaves.at(i) = ((double)numGt)/((double)fD.numNodes);
  }
  
  // subGraph only K steps out
  fD.subGraphKval = 4;
  getSubGraph(fD.numNodes,&fD.network,&fD.subGraphK,fD.subGraphKval);
  fD.subGraphKmax = 0;
  for(i=0; i<fD.numNodes; i++)
    if(fD.subGraphKmax < fD.subGraphK.at(i))
      fD.subGraphKmax = fD.subGraphK.at(i);
  
  // invGDistSD
  double mn=0,mnSq=0,d;
  tot=0;
  for(i=0; i<fD.numNodes; i++){
    for(j=(i+1); j<fD.numNodes; j++){
      d=1.0/(1.0+fD.gDist.at(i*fD.numNodes + j));
      mn+=d;
      mnSq+=d*d;
      tot++;
    }
  }
  mn/=(double)(tot);
  mnSq/=(double)(tot);
  fD.invGDistSD = std::sqrt(((double)(tot/(tot-1)))*(mnSq-mn*mn));
  
  // expInvGDistSD
  fD.expInvGDistSD.clear();
  fD.expInvGDistSD.reserve(fD.numNodes*fD.numNodes);
  for(i=0; i<fD.numNodes; i++){
    for(j=0; j<fD.numNodes; j++){
      d=fD.gDist.at(i*fD.numNodes+j);
      fD.expInvGDistSD.push_back(std::exp((1.0/(1.0+d))/fD.invGDistSD));
    }
  }

  // gDistSD
  tot = 0;
  for(i = 0; i < fD.numNodes; ++i){
    for(j = (i+1); j < fD.numNodes; ++j){
      d = fD.gDist.at(i*fD.numNodes + j);
      mn += d;
      mnSq += d*d;
      ++tot;
    }
  }
  mn /= (double)(tot);
  mnSq /= (double)(tot);
  fD.gDistSD = std::sqrt(((double)(tot/(tot-1)))*(mnSq-mn*mn));

  // expGDistSD
  fD.expGDistSD.clear();
  fD.expGDistSD.reserve(fD.numNodes*fD.numNodes);
  for(i=0; i<fD.numNodes; i++){
    for(j=0; j<fD.numNodes; j++){
      d=fD.gDist.at(i*fD.numNodes+j);
      fD.expGDistSD.push_back(std::exp(-d*d/(2*fD.gDistSD*fD.gDistSD)));
    }
  }

  // logGDist
  fD.logGDist.clear();
  fD.logGDist.reserve(fD.numNodes*fD.numNodes);
  for(i=0; i<fD.numNodes; i++){
    for(j=0; j<fD.numNodes; j++){
      d=fD.gDist.at(i*fD.numNodes+j);
      fD.logGDist.push_back(std::log(2.0+d));
    }
  }


  // half plane data depth
  fD.hpdd.clear();
  fD.hpdd.reserve(fD.numNodes);
  for(i = 0; i < fD.numNodes; ++i){
    fD.hpdd.push_back(halfPlaneDepth(fD.centroidsLong.at(i),
				     fD.centroidsLat.at(i),
				     fD.numNodes,
				     fD.centroidsLong,
				     fD.centroidsLat));
  }
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::nextPoint(){
  modelGen.modFill(sD,tD,fD,dD);
  modelGen.infProbs(sD,tD,fD,dD);

  // std::vector<double> infProbs = modelGen.infProbs();

  // double n = sD.numNotInfec;
  // double mn = std::accumulate(infProbs.begin(),infProbs.end(),0.0);
  // mn/=n;
  // double sd = std::accumulate(infProbs.begin(),infProbs.end(),0.0,
  // 			      [&mn](const double a, const double b){
  // 				return a + (b-mn)*(b-mn);
  // 			      });
  // sd/=n;
  
  // std::stringstream ss;
  // ss << " t: " << sD.time << std::endl
  //    << " o: " << value() << std::endl
  //    << " n: " << n << std::endl
  //    << "nn: " << sD.numNotInfec << std::endl
  //    << "ei: " << mn*n << std::endl
  //    << "mn: " << mn << std::endl
  //    << "sd: " << sd << std::endl;

  // njm::message(ss.str());

  // std::vector<double> vals = modelGen.probs;

  // mn = std::accumulate(vals.begin(),vals.end(),0.0);

  // njm::message("raw probs: " + njm::toString(mn,""));
  

					       
  nextPoint(modelGen.infProbs());
}


template<class MG,
	 class ME>
void System<MG,
	    ME>
::nextPoint(const std::vector<double> & infProbs){
  int i;
  int node,numNewInf=0;
  sD.newInfec.clear();
  double R,Rtot = 0;
  double propTot = 0;
  for(i=0; i<sD.numNotInfec; i++){
    R = njm::runif01();
    Rtot += R;
    propTot += infProbs.at(i);
    if(R < infProbs.at(i)){
      node = sD.notInfec.at(i);
      sD.infected.push_back(node);
      sD.newInfec.push_back(node);
      sD.notInfec.at(i)=fD.numNodes; // assign it the max value
      numNewInf++;
    }
  }


  // std::stringstream ss;
  // ss << "propTot: " << propTot << std::endl
  //    << "   Rtot: " << Rtot << std::endl
  //    << "  Rmean: " << Rtot/double(sD.numNotInfec) << std::endl
  //    << " numNew: " << numNewInf << std::endl;
  // njm::message(ss.str());
  
  std::sort(sD.infected.begin(),sD.infected.end());
  std::sort(sD.notInfec.begin(),sD.notInfec.end());
  sD.numInfected+=numNewInf;
  sD.numNotInfec-=numNewInf;

  // update time infected
  for(i=0; i<sD.numInfected; i++)
    ++sD.timeInf.at(sD.infected.at(i));

  // remove nodes that became infected
  sD.notInfec.erase(sD.notInfec.begin() + sD.numNotInfec, sD.notInfec.end());

  // update history
  sD.history.push_back(sD.status);
  // wipe treatment vectors
  tD.pPast = tD.p;
  tD.aPast = tD.a;
  std::fill(tD.p.begin(),tD.p.end(),0);
  std::fill(tD.a.begin(),tD.a.end(),0);

  updateStatus();
  
  sD.time++; // turn the calendar
}



template<class MG,
	 class ME>
void System<MG,
	    ME>::updateStatus(){
  int i,j,k,isInf;
  for(i=0,j=0,k=0; i<fD.numNodes; i++){
    if(j == sD.numInfected)
      isInf=0;
    else if(k == sD.numNotInfec)
      isInf = 1;
    else if(sD.infected.at(j) < sD.notInfec.at(k))
      isInf = 1;
    else
      isInf = 0;
    
    if(isInf){
      j++;
      if(tD.a.at(i))
	sD.status.at(i) = 3;
      else
	sD.status.at(i) = 2;
    }
    else{
      k++;
      if(tD.p.at(i))
	sD.status.at(i) = 1;
      else
	sD.status.at(i) = 0;
    }
  }
}


template<class MG,
	 class ME>
double System<MG,
	      ME>::value(){
  return ((double)sD.numInfected)/((double)fD.numNodes);
}
