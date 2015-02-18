#include "system.hpp"

template class System<GravityModel,GravityParam,
		      GravityModel,GravityParam>;

template class System<GravityModel,GravityParam,
		      RangeModel,RangeParam>;

template class System<GravityModel,GravityParam,
		      CaveModel,CaveParam>;

template class System<RangeModel,RangeParam,
		      RangeModel,RangeParam>;

template class System<CaveModel,CaveParam,
		      CaveModel,CaveParam>;

template class System<EbolaModel,EbolaParam,
		      EbolaModel,EbolaParam>;


template <class MG, class MPG,
	  class ME, class MPE>
System<MG, MPG,
       ME, MPE>::System(){
  initialize();
}


template <class MG, class MPG,
	  class ME, class MPE>
System<MG, MPG,
       ME, MPE>::System(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD,
			const MG & modelGen,
			const ME & modelEst,
			const MPG & paramGen,
			const MPE & paramEst){
  this->sD_r = sD;
  this->tD_r = tD;
  this->fD = fD;
  this->dD_r = dD;
  this->modelGen = modelGen;
  this->modelEst = modelEst;
  this->paramGen_r = paramGen;
  this->paramEst_r = paramEst;
  reset();
}


template <class MG, class MPG,
	  class ME, class MPE>
System<MG, MPG,
       ME, MPE>::System(const std::string file){
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

  reset();
}


template <class MG, class MPG,
	  class ME, class MPE>
void System<MG, MPG,
	    ME, MPE>::reset(){
  sD = sD_r;
  tD = tD_r;
  dD = dD_r;

  paramGen = paramGen_r;
  paramEst = paramEst_r;
}


template <class MG, class MPG,
	  class ME, class MPE>
void System<MG, MPG,
	    ME, MPE>::checkPoint(){
  sD_r = sD;
  tD_r = tD;
  dD_r = dD;
  
  paramGen_r = paramGen;
  paramEst_r = paramEst;
}


template <class MG, class MPG,
	  class ME, class MPE>
void System<MG, MPG,
	    ME, MPE>::initialize(){
  njm::fromFile(fD.fips,njm::sett.srcExt("fips.txt"));
  fD.numNodes = fD.fips.size();
  njm::fromFile(fD.dist,njm::sett.srcExt("d.txt"));
  njm::fromFile(fD.caves,njm::sett.srcExt("caves.txt"));
  njm::fromFile(fD.covar,njm::sett.srcExt("xcov.txt"));
  fD.numCovar = ((int)fD.covar.size())/fD.numNodes;
  njm::fromFile(fD.network,njm::sett.srcExt("network.txt"));
  njm::fromFile(fD.centroidsLong,njm::sett.srcExt("centroidsLong.txt"));
  njm::fromFile(fD.centroidsLat,njm::sett.srcExt("centroidsLat.txt"));
  njm::fromFile(fD.subGraph,njm::sett.srcExt("subGraph.txt"));
  njm::fromFile(fD.betweenness,njm::sett.srcExt("betweenness.txt"));

  std::vector<double> start;
  njm::fromFile(start,njm::sett.srcExt("startingLocations.txt"));

  
  sD_r.time=0;
  sD_r.history.clear();
  sD_r.status.resize(fD.numNodes);
  std::fill(sD_r.status.begin(),sD_r.status.end(),0);
  int i;
  for(i=0; i<(int)start.size(); i++)
    sD_r.status.at(start.at(i))=2;
  sD_r.numInfected = 0;
  sD_r.numNotInfec = 0;    
  for(i=0; i<fD.numNodes; i++)
    if(sD_r.status.at(i)==2){
      sD_r.infected.push_back(i);
      sD_r.timeInf.push_back(1);
      sD_r.numInfected++;
    }
    else{
      sD_r.notInfec.push_back(i);
      sD_r.timeInf.push_back(0);      
      sD_r.numNotInfec++;
    }
  sD_r.newInfec = sD_r.infected;
  tD_r.a.resize(fD.numNodes);
  std::fill(tD_r.a.begin(),tD_r.a.end(),0);
  tD_r.p.resize(fD.numNodes);
  std::fill(tD_r.p.begin(),tD_r.p.end(),0);
  tD_r.aPast.resize(fD.numNodes);
  std::fill(tD_r.aPast.begin(),tD_r.aPast.end(),0);
  tD_r.pPast.resize(fD.numNodes);
  std::fill(tD_r.pPast.begin(),tD_r.pPast.end(),0);


  // only start treatment at time trtStart and on  
  njm::fromFile(fD.trtStart,njm::sett.srcExt("trtStart.txt"));
  // only update every period steps
  njm::fromFile(fD.period,njm::sett.srcExt("period.txt"));
  // final time step in simulation
  njm::fromFile(fD.finalT,njm::sett.srcExt("finalT.txt"));

  preCompData();

  paramGen_r.load();
  modelGen.load(sD_r,tD_r,fD,dD_r,paramGen_r);

  reset();
}


template <class MG, class MPG,
	  class ME, class MPE>
void System<MG,MPG,
	    ME,MPE>::preCompData(){
  int i,j,tot;
  
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


template <class MG, class MPG,
	  class ME, class MPE>
void System<MG, MPG,
	    ME, MPE>::nextPoint(){
  modelGen.infProbs(sD,tD,fD,dD,paramGen);
  nextPoint(paramGen.infProbs);
}


template<class MG, class MPG,
	 class ME, class MPE>
void System<MG, MPG,
	    ME, MPE>
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



template<class MG, class MPG,
	 class ME, class MPE>
void System<MG, MPG,
	    ME, MPE>::updateStatus(){
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


template<class MG, class MPG,
	 class ME, class MPE>
double System<MG, MPG,
	      ME, MPE>::value(){
  return ((double)sD.numInfected)/((double)fD.numNodes);
}



////////////////////////////////////////////////////////////////////////////////
// System Light



template class SystemLight<GravityModel,GravityParam>;


template<class M, class MP>
SystemLight<M, MP>::SystemLight(const SimData & sD,
				const TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD,
				const M & modelGen,
				const MP & paramGen){
  this->sD_r = sD;
  this->tD_r = tD;
  this->fD = fD;
  this->dD_r = dD;
  this->modelGen = modelGen;
  this->paramGen_r = paramGen;

  reset();
}


template<class M, class MP>
void SystemLight<M, MP>::reset(){
  sD = sD_r;
  tD = tD_r;
  dD = dD_r;

  paramGen = paramGen_r;
}



template<class M, class MP>
void
SystemLight<M, MP>::nextPoint(const int isFinal){
  int i;
  for(i=0; i<sD.numInfected; i++)
    sD.timeInf.at(sD.infected.at(i))++;
  
  int node,numNewInf=0;
  sD.newInfec.clear();
  for(i=0; i<sD.numNotInfec; i++){
    if(njm::runif01() < paramGen.infProbs.at(i)){
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
  
  sD.time++; // turn the calendar

  if(!isFinal)
    modelGen.update(sD,tD,fD,dD,paramGen);
}



template<class M, class MP>
double SystemLight<M, MP>::value(){
  return ((double)sD.numInfected)/((double)fD.numNodes);
}


