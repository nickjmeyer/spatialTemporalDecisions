#include "bayesP.hpp"

std::vector<double> getStats(const std::vector<std::vector<int> > & h,
			     const SimData & sD,
			     const FixedData & fD){
  std::vector<double> stats;
  // total infected
  stats.push_back(sD.numInfected);

  int i,j,sum = 0;
  for(i = 1; i < (int)h.size(); ++i){
    sum = 0;
    for(j = 0; j < fD.numNodes; ++j)
      if(h.at(i).at(j) >= 2 && h.at(i-1).at(j) < 2)
  	++sum;

    stats.push_back(sum);
  }


  // average year of infection
  sum = 0;
  for(i = 0; i < (int)h.size(); ++i){
    for(j = 0; j < fD.numNodes; ++j){
      if((i == 0 && h.at(i).at(j) >= 2) ||
	 (i > 0 && h.at(i).at(j) >= 2 && h.at(i-1).at(j) < 2))
	sum += i;
    }
  }
  stats.push_back(((double)sum)/((double)sD.numInfected));


  // average lat long for infected
  double cLong,cLat;
  cLong = cLat = 0;
  for(i = 0 ; i < sD.numInfected; ++i){
    cLong += fD.centroidsLong.at(sD.infected.at(i));
    cLat += fD.centroidsLat.at(sD.infected.at(i));
  }
  stats.push_back(cLong/((double)sD.numInfected));
  stats.push_back(cLat/((double)sD.numInfected));

  // average spread distance from start for infected
  double sDist;
  int start = -1,numStart = 0;
  for(i = 0; i < fD.numNodes; ++i){
    if(h.at(0).at(i) >= 2){
      if(numStart == 0)
	start = i;
      ++numStart;
    }
  }
  if(numStart > 1)
    std::cout << "warning...multiple starting locations..."
	      << std::endl;
  sDist = 0;
  for(i = 0; i < sD.numInfected; ++i){
    sDist += fD.dist.at(sD.infected.at(i)*fD.numNodes + start);
  }
  sDist/=(double)sD.numInfected;
  stats.push_back(sDist);

  // min lat long for infected
  cLong = std::numeric_limits<double>::max();
  cLat = std::numeric_limits<double>::max();
  for(i = 0; i < sD.numInfected; ++i){
    if(fD.centroidsLong.at(sD.infected.at(i)) < cLong)
      cLong = fD.centroidsLong.at(sD.infected.at(i));
    if(fD.centroidsLat.at(sD.infected.at(i)) < cLat)
      cLat = fD.centroidsLat.at(sD.infected.at(i));
  }
  stats.push_back(cLong);
  stats.push_back(cLat);


  // max lat long for infected
  cLong = std::numeric_limits<double>::lowest();
  cLat = std::numeric_limits<double>::lowest();
  for(i = 0; i < sD.numInfected; ++i){
    if(fD.centroidsLong.at(sD.infected.at(i)) > cLong)
      cLong = fD.centroidsLong.at(sD.infected.at(i));
    if(fD.centroidsLat.at(sD.infected.at(i)) > cLat)
      cLat = fD.centroidsLat.at(sD.infected.at(i));
  }
  stats.push_back(cLong);
  stats.push_back(cLat);

  
  // max dist from starting location
  sDist = std::numeric_limits<double>::lowest();
  for(i = 0; i < sD.numInfected; ++i)
    if(fD.dist.at(sD.infected.at(i)*fD.numNodes + start) > sDist)
      sDist = fD.dist.at(sD.infected.at(i)*fD.numNodes + start);
  stats.push_back(sDist);


  return stats;
}



int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfModel GM;
  typedef GravityTimeInfParam GP;
  typedef GM EM;
  typedef GP EP;

  typedef System<GM,GP,EM,EP> S;


  S s1Obs("obsData.txt");

  std::vector<std::vector<int> > h;
  h = s1Obs.sD.history;
  h.push_back(s1Obs.sD.status);

  std::vector<std::string> names = {"n_inf","n_inf_2007","n_inf_2008",
  				    "n_inf_2009","n_inf_2010","n_inf_2011",
  				    "n_inf_2012","n_inf_2013","mean_year",
  				    "mean_long","mean_lat",
  				    "mean_dist_from_start",       
  				    "min_long","min_lat",       
  				    "max_long","max_lat",
  				    "max_dist_from_start"};

  njm::toFile(names,njm::sett.datExt("obsStats_",".txt"),
  	      std::ios_base::out);
  njm::toFile(getStats(h,s1Obs.sD,s1Obs.fD),
	      njm::sett.datExt("obsStats_",".txt"));

  s1Obs.modelGen.mcmc.load(s1Obs.sD.history,s1Obs.sD.status,s1Obs.fD);
  s1Obs.modelGen.mcmc.sample(20000,10000);

  std::vector< std::vector<double> > stats;
  
  S s1;
  s1.modelGen = s1Obs.modelGen;
  int r,t,R,T;
  R = 10000;
  T = s1Obs.sD.time;
  for(r = 0; r < R; ++r){
    s1.modelGen.mcmc.samples.setRand();
    s1.paramGen_r.putPar(s1.modelGen.mcmc.samples.getPar());
    s1.paramEst_r.putPar(s1.paramGen_r.getPar());
    
    s1.reset();
    
    for(t = 0; t < T; ++t)
      s1.nextPoint();
    h = s1.sD.history;
    h.push_back(s1.sD.status);

    stats.push_back(getStats(h,s1.sD,s1.fD));
  }

  njm::toFile(names,njm::sett.datExt("sampStats_",".txt"),
  	      std::ios_base::out);
  njm::toFile(njm::toString(stats,"\n",""),
	      njm::sett.datExt("sampStats_timeInf_",".txt"));


  
  stats.clear();

  
  typedef GravityModel GM2;
  typedef GravityParam GP2;
  typedef GM2 EM2;
  typedef GP2 EP2;


  System<GM2,GP2,EM2,EP2> s2Obs("obsData.txt");
  
  s2Obs.modelGen.mcmc.load(s2Obs.sD.history,s2Obs.sD.status,s2Obs.fD);
  s2Obs.modelGen.mcmc.sample(20000,10000);

  System<GM2,GP2,EM2,EP2> s2;

  s2.modelGen = s2Obs.modelGen;

  for(r = 0; r < R; ++r){
    s2.modelGen.mcmc.samples.setRand();
    s2.paramGen_r.putPar(s2.modelGen.mcmc.samples.getPar());
    s2.paramEst_r.putPar(s2.paramGen_r.getPar());
    
    s2.reset();
    
    for(t = 0; t < T; ++t)
      s2.nextPoint();
    h = s2.sD.history;
    h.push_back(s2.sD.status);

    stats.push_back(getStats(h,s2.sD,s2.fD));
  }

  njm::toFile(names,njm::sett.datExt("sampStats_",".txt"),
  	      std::ios_base::out);
  njm::toFile(njm::toString(stats,"\n",""),
	      njm::sett.datExt("sampStats_",".txt"));
  

  // njm::sett.clean();
  return 0;
}
