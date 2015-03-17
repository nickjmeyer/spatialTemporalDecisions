#include "tuneTrtProp.hpp"


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

  std::vector<double> agent = {0.0,1.0,2.0,3.0};
  std::vector<double> propTrt = {0.01,0.03,0.06,0.1};
  

  FFX ffx;

  ffx.addFactor("agent",agent);
  ffx.addFactor("propTrt",propTrt);

  ffx.addStat("value");
  ffx.addStat("time");

  ffx.setReps(4);

  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef MG ME;
  typedef PG PE;

  typedef System<MG,PG,ME,PE> S;

  typedef NoTrt<ME,PE> NT;
  typedef ProximalAgent<ME,PE> PA;
  typedef MyopicAgent<ME,PE> MA;
  
  typedef ToyFeatures2<ME,PE> F;
  typedef RankAgent<F,ME,PE> RA;

  typedef M1SpOptim<S,RA,ME,PE> SPO;

  typedef VanillaRunnerNS<S,NT> R_NT;
  typedef VanillaRunnerNS<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunnerNS<S,RA,SPO> R_RA;


  S s;
  s.modelGen.fitType = MLE;  // for speed
  s.modelEst.fitType = MLE;  // for speed

  NT nt;
  PA pa;
  MA ma;
  RA ra;
  ra.reset();

  SPO spo;

  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_RA r_ra;
  

  
  // this is an experiment to set up the tuning, so no tuning necessary
  spo.tp.tune = 0;

  double value;
  double agentInd;
  int done = 0;
  int i, M = ffx.maxInd;
  int tick,tock;
  for(i = 0; i < M; ++i){
    s.fD.propTrt = ffx.getSett("propTrt",i);
    agentInd = ffx.getSett("agent",i);

    tick = std::time(NULL);
    if(agentInd < 0.1)
      value = r_nt.run(s,nt,300,s.fD.finalT);
    else if(0.9 < agentInd && agentInd < 1.1)
      value = r_pa.run(s,pa,300,s.fD.finalT);
    else if(1.9 < agentInd && agentInd < 2.1)
      value = r_ma.run(s,ma,300,s.fD.finalT);
    else if(2.9 < agentInd && agentInd < 3.1)
      value = r_ra.run(s,ra,spo,300,s.fD.finalT);
    else
      throw(1);
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

