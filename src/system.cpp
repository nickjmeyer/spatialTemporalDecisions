#include "system.hpp"

template class System<GravityTimeInfExpCavesModel,
		      GravityTimeInfExpCavesModel>;

template class System<GravityTimeInfExpCavesModel,
		      GravityTimeInfExpModel>;

template class System<GravityTimeInfExpCavesModel,
		      GravityTimeInfModel>;

template class System<GravityTimeInfExpCavesModel,
		      GravityModel>;

template class System<GravityTimeInfExpCavesModel,
		      RangeModel>;

template class System<GravityTimeInfExpCavesModel,
		      RadiusModel>;

template class System<GravityTimeInfExpCavesModel,
		      CaveModel>;

template class System<GravityTimeInfExpCavesModel,
		      MultiModel>;



template class System<GravityTimeInfExpModel,
		      GravityTimeInfExpModel>;

template class System<GravityTimeInfModel,
		      GravityTimeInfModel>;

template class System<GravityModel,
		      GravityModel>;

template class System<RangeModel,
		      RangeModel>;

template class System<RadiusModel,
		      RadiusModel>;

template class System<CaveModel,
		      CaveModel>;

template class System<MultiModel,
		      MultiModel>;


template <class MG,
	  class ME>
System<MG,
       ME>::System(){
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
  for(j = 0; j < fD.numNodes; j++){
    tD_r.a.at(j) = tD_r.p.at(j) = 0;
    if(sD_r.status.at(j) == 1)
      tD_r.p.at(j) = 1;
    else if(sD_r.status.at(j) == 3)
      tD_r.a.at(j) = 1;
  }

  // past treatments
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

  revert();
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::reset(const std::vector<int> & ind){
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
  modelGen.load(sD_r,tD_r,fD,dD_r);

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
  njm::fromFile(fD.dist,njm::sett.srcExt("d.txt"));
  njm::fromFile(fD.caves,njm::sett.srcExt("caves.txt"));
  njm::fromFile(fD.covar,njm::sett.srcExt("xcov.txt"));
  fD.numCovar = ((int)fD.covar.size())/fD.numNodes;
  njm::fromFile(fD.network,njm::sett.srcExt("network.txt"));
  
  njm::fromFile(fD.centroidsMdsLong,njm::sett.srcExt("centroidsLong.txt"));
  njm::fromFile(fD.centroidsMdsLat,njm::sett.srcExt("centroidsLat.txt"));
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

  modelGen.setType(MCMC);
  modelEst.setType(MCMC);

  modelGen_r.getPar()->load();
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::preCompData(){
  int i,j,tot;

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
  
  // invDistSD
  double mn=0,mnSq=0,d;
  tot=0;
  for(i=0; i<fD.numNodes; i++){
    for(j=(i+1); j<fD.numNodes; j++){
      d=1.0/(1.0+fD.dist.at(i*fD.numNodes + j));
      mn+=d;
      mnSq+=d*d;
      tot++;
    }
  }
  mn/=(double)(tot);
  mnSq/=(double)(tot);
  fD.invDistSD = std::sqrt(((double)(tot/(tot-1)))*(mnSq-mn*mn));
  
  // expInvDistSD
  fD.expInvDistSD.clear();
  fD.expInvDistSD.reserve(fD.numNodes*fD.numNodes);
  for(i=0; i<fD.numNodes; i++){
    for(j=0; j<fD.numNodes; j++){
      d=fD.dist.at(i*fD.numNodes+j);
      fD.expInvDistSD.push_back(std::exp((1.0/(1.0+d))/fD.invDistSD));
    }
  }

  // logDist
  fD.logDist.clear();
  fD.logDist.reserve(fD.numNodes*fD.numNodes);
  for(i=0; i<fD.numNodes; i++){
    for(j=0; j<fD.numNodes; j++){
      d=fD.dist.at(i*fD.numNodes+j);
      fD.logDist.push_back(std::log(2.0+d));
    }
  }
}


template <class MG,
	  class ME>
void System<MG,
	    ME>::nextPoint(){
  modelGen.infProbs(sD,tD,fD,dD);
  nextPoint(modelGen.getPar()->getInfProbs());
}


template<class MG,
	 class ME>
void System<MG,
	    ME>
::nextPoint(const std::vector<double> & infProbs){
  int i;
  for(i=0; i<sD.numInfected; i++)
    sD.timeInf.at(sD.infected.at(i))++;
  
  int node,numNewInf=0;
  sD.newInfec.clear();
  for(i=0; i<sD.numNotInfec; i++){
    if(njm::runif01() < infProbs.at(i)){
      node = sD.notInfec.at(i);
      sD.infected.push_back(node);
      sD.newInfec.push_back(node);
      sD.notInfec.at(i)=fD.numNodes; // assign it the max value
      numNewInf++;
    }
  }
  std::sort(sD.infected.begin(),sD.infected.end());
  std::sort(sD.notInfec.begin(),sD.notInfec.end());
  sD.numInfected+=numNewInf;
  sD.numNotInfec-=numNewInf;

  sD.notInfec.erase(sD.notInfec.begin() + sD.numNotInfec, sD.notInfec.end());

  sD.history.push_back(sD.status);
  updateStatus();
  
  // wipe treatment vectors
  tD.pPast = tD.p;
  tD.aPast = tD.a;
  std::fill(tD.p.begin(),tD.p.end(),0);
  std::fill(tD.a.begin(),tD.a.end(),0);

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



