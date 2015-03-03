#include "tuneSp.hpp"


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

  std::vector<double> Avals = {30,50};
  std::vector<double> Bvals = {0.3,0.5};
  std::vector<double> Cvals = {2.0,10.0};
  std::vector<double> Tvals = {0.2,1.0,2.0};
  std::vector<double> Lvals = {0.5,1.0,2.0};

  FFX ffx;

  ffx.addFactor("A",Avals);
  ffx.addFactor("B",Bvals);
  ffx.addFactor("C",Cvals);
  ffx.addFactor("T",Tvals);
  ffx.addFactor("L",Lvals);

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
  AR ar;
  SPO spo;
  SPR spr;

  FR fr;

  ar.reset();
  njm::message("Fit only: " + njm::toString(fr.run(s,ar,300,s.fD.finalT)));

  // for speed
  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  // this is an experiment to set up the tuning, so no tuning necessary
  spo.tp.tune = 0;

  double value;
  int done = 0;
  int i, M = ffx.maxInd;
  int tick,tock;
  for(i = 0; i < M; ++i){
    spo.tp.A = ffx.getSett("A",i);
    spo.tp.B = ffx.getSett("B",i);
    spo.tp.C = ffx.getSett("C",i);
    spo.tp.t = ffx.getSett("T",i);
    spo.tp.ell = ffx.getSett("L",i);

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

