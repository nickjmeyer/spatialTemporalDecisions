#include "system.hpp"

template class System<GravityModel,GravityParam>;

template class System<EbolaModel,EbolaParam>;


template<class Model, class ModelParam>
System<Model, ModelParam>::System(){
  initialize();
}


template<class Model, class ModelParam>
System<Model, ModelParam>::System(const SimData & sD,
				  const TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  const Model & model,
				  const ModelParam & genParam,
				  const ModelParam & estParam){
  this->sD_r = sD;
  this->tD_r = tD;
  this->fD = fD;
  this->dD_r = dD;
  this->model = model;
  this->genParam_r = genParam;
  this->estParam_r = estParam;

  reset();
}


template<class Model, class ModelParam>
void System<Model, ModelParam>::reset(){
  sD = sD_r;
  tD = tD_r;
  dD = dD_r;

  genParam = genParam_r;
  estParam = estParam_r;
}


template<class Model, class ModelParam>
void System<Model, ModelParam>::checkPoint(){
  sD_r = sD;
  tD_r = tD;
  dD_r = dD;
  
  genParam_r = genParam;
  estParam_r = estParam;
}


template<class Model, class ModelParam>
void System<Model, ModelParam>::initialize(){
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

  genParam_r.load();
  model.load(sD_r,tD_r,fD,dD_r,genParam_r);

  reset();

  preCompData();

}


template<class Model, class ModelParam>
void System<Model,ModelParam>::preCompData(){
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


template<class Model, class ModelParam>
void System<Model, ModelParam>::nextPoint(){
  model.infProbs(sD,tD,fD,dD,genParam);
  nextPoint(genParam.infProbs);
}


template<class Model, class ModelParam>
void System<Model, ModelParam>::nextPoint(const std::vector<double> & infProbs){
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



template<class Model, class ModelParam>
void System<Model,ModelParam>::updateStatus(){
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


template<class Model, class ModelParam>
double System<Model, ModelParam>::value(){
  return ((double)sD.numInfected)/((double)fD.numNodes);
}







template class SystemLight<GravityModel,GravityParam>;


template<class Model, class ModelParam>
SystemLight<Model, ModelParam>::SystemLight(const SimData & sD,
					    const TrtData & tD,
					    const FixedData & fD,
					    const DynamicData & dD,
					    const Model & model,
					    const ModelParam & genParam){
  this->sD_r = sD;
  this->tD_r = tD;
  this->fD = fD;
  this->dD_r = dD;
  this->model = model;
  this->genParam_r = genParam;

  reset();
}


template<class Model, class ModelParam>
void SystemLight<Model, ModelParam>::reset(){
  sD = sD_r;
  tD = tD_r;
  dD = dD_r;

  genParam = genParam_r;
}



template<class Model, class ModelParam>
void
SystemLight<Model, ModelParam>::nextPoint(const int isFinal){
  int i;
  for(i=0; i<sD.numInfected; i++)
    sD.timeInf.at(sD.infected.at(i))++;
  
  int node,numNewInf=0;
  sD.newInfec.clear();
  for(i=0; i<sD.numNotInfec; i++){
    if(njm::runif01() < genParam.infProbs.at(i)){
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
    model.update(sD,tD,fD,dD,genParam);
}



template<class Model, class ModelParam>
double SystemLight<Model, ModelParam>::value(){
  return ((double)sD.numInfected)/((double)fD.numNodes);
}


