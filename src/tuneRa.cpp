#include "tuneRa.hpp"


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

  std::vector<double> chunkVals = {1,2,3};
  std::vector<double> scaleVals = {1,2,4,8,10};

  FFX ffx;

  ffx.addFactor("chunk",chunkVals);
  ffx.addFactor("scale",scaleVals);

  ffx.addStat("value");
  ffx.addStat("time");

  ffx.setReps(8);

  typedef GravityTimeInfModel MG;
  typedef GravityTimeInfParam PG;
  
  typedef MG ME;
  typedef PG PE;

  typedef System<MG,PG,ME,PE> S;
  
  typedef ToyFeatures2<ME,PE> F;
  
  typedef RankAgent<F,ME,PE> AR;

  typedef M1SpOptim<S,AR,ME,PE> SPO;

  typedef OptimRunnerNS<S,AR,SPO> SPR;

  typedef FitOnlyRunner<S,AR> FR;

  S s;
  s.modelGen.fitType = MLE;  // for speed
  s.modelEst.fitType = MLE;  // for speed
  
  AR ar;
  ar.reset();

  SPO spo;
  
  SPR spr;

  FR fr;

  // baseline value
  njm::message(" Fit only: " + njm::toString(fr.run(s,ar,300,s.fD.finalT)));
  
  // this is an experiment done after tuning on grid 100, so we know
  // good values for the grid 100
  spo.tp.tune = 0;

  // original agent tp
  njm::message("Original: " + njm::toString(spr.run(s,ar,spo,300,s.fD.finalT)));

  double value;
  int done = 0;
  int i, M = ffx.maxInd;
  int tick,tock;
  for(i = 0; i < M; ++i){
    ar.tp.numChunks = ffx.getSett("chunk",i);
    ar.tp.jitterScale = ffx.getSett("scale",i);

    tick = std::time(NULL);
    value = spr.run(s,ar,spo,300,s.fD.finalT);
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

