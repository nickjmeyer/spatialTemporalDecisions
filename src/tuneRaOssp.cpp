#include "tuneRaOssp.hpp"


FFX::FFX(){
  numFactor = 0;
  numStat = 0;
  numCombo = 1;
  numReps = 1;
}


void FFX::setReps(const int num){
#pragma omp critical(FFX)
  {
    numReps = (num > 0 ? num : 1);
    maxInd = numReps * numCombo;
  }
}


void FFX::addFactor(const std::string & f,
		    const std::vector<double> & fVals){
#pragma omp critical(FFX)
  {
    factors.push_back(f);
    values.push_back(fVals);
    maxSett.push_back((int)fVals.size());

    numCombo *= (int)fVals.size();
    ++numFactor;

    maxInd = numReps * numCombo;
  }
}

void FFX::addStat(const std::string & s){
#pragma omp critical(FFX)
  {
    stats.push_back(s);
    ++numStat;
  }
}


std::vector<double> FFX::getSett(const int num) const {
  if(num >= (numCombo * numReps)){
    std::cout << "invalid setting index" << std::endl;
    throw(1);
  }

  int i, maxInd = numCombo, ind = num % numCombo;
  std::vector<double> sett;
  for(i = 0 ; i < numFactor; i++){
    maxInd /= (int)maxSett.at(i);
    sett.push_back(values.at(i).at(ind/maxInd));
    ind %= maxInd;
  }
  return sett;
}

double FFX::getSett(const std::string & f,const int num) const{
  std::vector<double> setting = getSett(num);
  int i;
  for(i = 0; i < numFactor; ++i)
    if(factors.at(i) == f)
      return setting.at(i);
  std::cout << "invalid factor name: " << f << std::endl;
  throw(1);
  return 0;
}


void FFX::addObs(const int num, const double & obs){
  std::vector<double> obsV;
  obsV.push_back(obs);
  addObs(num,obsV);
}

void FFX::addObs(const int num, const std::vector<double> & obs){
#pragma omp critical(FFX)
  {
    if(numStat != (int)obs.size()){
      std::cout << "invalid number of stats: " << (int)obs.size()
		<< std::endl;
      throw(1);
    }

    std::vector<double> insObs;
    insObs.push_back(num % numCombo); // insert combo number for convenience
    std::vector<double> sett = getSett(num);
    insObs.insert(insObs.end(),sett.begin(),sett.end());
    insObs.insert(insObs.end(),obs.begin(),obs.end());

    allObs.push_back(insObs);
  }
}


void FFX::saveObs(const std::string & file) const{
#pragma omp critical(FFX)
  {
    // clears file and add header
    std::vector<std::string> header = {"combo"};
    header.insert(header.end(),factors.begin(),factors.end());
    header.insert(header.end(),stats.begin(),stats.end());
  
    njm::toFile(njm::toString(header," ","\n"),file,std::ios_base::out);
    // append results
    njm::toFile(njm::toString(allObs,"\n",""),file,std::ios_base::app);
  }
}



int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  std::vector<double> scaleVals = {2,4,8};
  std::vector<double> nVals = {100,500,1000};

  FFX ffx;

  ffx.addFactor("scale",scaleVals);
  ffx.addFactor("N",nVals);

  ffx.addStat("value");
  ffx.addStat("time");

  ffx.setReps(4);

  typedef GravityTimeInfExpCavesModel MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef ToyFeatures2<ME> F;
  typedef OsspAgent<ME> OA;

  typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef OptimRunner<S,OA,OSSPO> R_OA;

  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  OA oa;

  OSSPO osspo;

  R_OA r_oa;  

  double value;
  int done = 0;
  int i, M = ffx.maxInd;
  int tick,tock;
  for(i = 0; i < M; ++i){
    osspo.tp.jitterScale = ffx.getSett("scale",i);
    osspo.tp.N = ffx.getSett("N",i);

    njm::randomSeed = unsigned(i / ffx.numCombo);

    tick = std::time(NULL);
    value = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
    tock = std::time(NULL);

    ffx.addObs(i,{value,((double)(tock-tick))/3600.0});
    ffx.saveObs(njm::sett.datExt("results_",".txt"));

    printf("\rFinished % 5d out of % 5d",++done,ffx.maxInd);
    fflush(stdout);
  }
  printf("\n");

  ffx.saveObs(njm::sett.datExt("results_",".txt"));

  return 0;
}

