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


template <class M, class P>
void runBayesP(const std::string & file, const int obs,
	       const int numSamples,const int numBurn,
	       const int numStats){
  typedef System<M,P,M,P> S;

  S sObs("obsData.txt");

  std::vector<std::vector<int> > h;
  h = sObs.sD.history;
  h.push_back(sObs.sD.status);

  std::vector<std::string> names = {"n_inf","n_inf_2007","n_inf_2008",
  				    "n_inf_2009","n_inf_2010","n_inf_2011",
  				    "n_inf_2012","n_inf_2013","mean_year",
  				    "mean_long","mean_lat",
  				    "mean_dist_from_start",       
  				    "min_long","min_lat",       
  				    "max_long","max_lat",
  				    "max_dist_from_start"};

  if(obs){
    njm::toFile(names,njm::sett.datExt("obsStats_",".txt"),
		std::ios_base::out);
    njm::toFile(getStats(h,sObs.sD,sObs.fD),
		njm::sett.datExt("obsStats_",".txt"));
  }

  sObs.modelGen.mcmc.load(sObs.sD.history,sObs.sD.status,sObs.fD);
  sObs.modelGen.mcmc.sample(numSamples,numBurn);

  sObs.modelGen.mcmc.samples.setMean();
  sObs.paramGen.putPar(sObs.modelGen.mcmc.samples.getPar());
  sObs.paramGen.save();

  std::vector< std::vector<double> > stats;
  
  S s;
  s.modelGen = sObs.modelGen;
  int r,t,R,T;
  R = numStats;
  T = sObs.sD.time;
  for(r = 0; r < R; ++r){
    s.modelGen.mcmc.samples.setRand();
    s.paramGen_r.putPar(s.modelGen.mcmc.samples.getPar());
    s.paramEst_r.putPar(s.paramGen_r.getPar());
    
    s.reset();
    
    for(t = 0; t < T; ++t)
      s.nextPoint();
    h = s.sD.history;
    h.push_back(s.sD.status);

    stats.push_back(getStats(h,s.sD,s.fD));
  }

  njm::toFile(names,njm::sett.datExt(file+"_",".txt"),
  	      std::ios_base::out);
  njm::toFile(njm::toString(stats,"\n",""),
	      njm::sett.datExt(file+"_",".txt"));

  std::vector<std::vector<double> > parSamp;
  int i;
  for(i = 0; i < (numSamples-numBurn); ++i){
    s.modelGen.mcmc.samples.setPar(i);
    parSamp.push_back(s.modelGen.mcmc.samples.getPar());
  }

  njm::toFile(njm::toString(parSamp,"\n",""),
	      njm::sett.datExt(file+"_param_",".txt"),
	      std::ios_base::out);
}



int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  int numSamples = 20000, numBurn = 10000, numStats = 10000;
  // int numSamples = 100, numBurn = 50, numStats = 50;

#pragma omp parallel sections			\
  shared(numSamples,numBurn,numStats)
  {
#pragma omp section
    {
      runBayesP<GravityModel,
		GravityParam>("sampStats_gravity",1,
			      numSamples,numBurn,numStats);
    }

#pragma omp section
    {
      runBayesP<GravityTimeInfModel,
		GravityTimeInfParam>("sampStats_timeInf",0,
				     numSamples,numBurn,numStats);
    }

// #pragma omp section
//     {
//       runBayesP<GravityTimeInfSqModel,
// 		GravityTimeInfSqParam>("sampStats_timeInfSq",0,
// 				       numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<GravityTimeInfSqrtModel,
// 		GravityTimeInfSqrtParam>("sampStats_timeInfSqrt",0,
// 					 numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<GravityTimeInfLogModel,
// 		GravityTimeInfLogParam>("sampStats_timeInfLog",0,
// 					numSamples,numBurn,numStats);
//     }

#pragma omp section
    {
      runBayesP<GravityTimeInfExpModel,
		GravityTimeInfExpParam>("sampStats_timeInfExp",0,
					numSamples,numBurn,numStats);
    }

#pragma omp section
    {
      runBayesP<GravityTimeInfExpCavesModel,
		GravityTimeInfExpCavesParam>("sampStats_timeInfExpCaves",0,
					     numSamples,numBurn,numStats);
    }

#pragma omp section
    {
      runBayesP<GravityTimeInfExpLCavesModel,
		GravityTimeInfExpLCavesParam>("sampStats_timeInfExpLCaves",0,
					      numSamples,numBurn,numStats);
    }

#pragma omp section
    {
      runBayesP<GravityTimeInfExpRCavesModel,
		GravityTimeInfExpRCavesParam>("sampStats_timeInfExpRCaves",0,
					      numSamples,numBurn,numStats);
    }
  
  }
  
  // njm::sett.clean();
  return 0;
}
