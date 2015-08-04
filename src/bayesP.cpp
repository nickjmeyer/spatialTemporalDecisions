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
    sDist += fD.gDist.at(sD.infected.at(i)*fD.numNodes + start);
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
    if(fD.gDist.at(sD.infected.at(i)*fD.numNodes + start) > sDist)
      sDist = fD.gDist.at(sD.infected.at(i)*fD.numNodes + start);
  stats.push_back(sDist);


  return stats;
}


template <class M>
void runBayesP(const std::string & file, const int obs,
	       const int numSamples,const int numBurn,
	       const int numStats){
  njm::resetSeed(0);

  typedef System<M,M> S;

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

  sObs.modelGen_r.mcmc.load(sObs.sD.history,sObs.sD.status,sObs.fD);
  sObs.modelGen_r.mcmc.sample(numSamples,numBurn,true);

  { // mean
    sObs.modelGen_r.mcmc.samples.setMean();
    std::vector<double> par = sObs.modelGen_r.mcmc.samples.getPar();
    sObs.modelGen_r.putPar(par.begin());
    sObs.modelGen_r.save();

    std::vector< std::vector<double> > stats;

    S s;
    Starts starts("startingLocations.txt");
    s.modelGen_r = sObs.modelGen_r;
    s.modelEst_r = s.modelGen_r;
    int r,t,R,T;
    R = numStats;
    T = sObs.sD.time;
    for(r = 0; r < R; ++r){
      s.modelGen_r.mcmc.samples.setRand();
      par = s.modelGen_r.mcmc.samples.getPar();
      s.modelGen_r.putPar(par.begin());
      s.modelEst_r.putPar(par.begin());

      s.reset(starts[r]);

      for(t = 0; t < T; ++t)
	s.nextPoint();
      h = s.sD.history;
      h.push_back(s.sD.status);

      stats.push_back(getStats(h,s.sD,s.fD));
    }

    njm::toFile(names,njm::sett.datExt("sampStats_mean_"+file+"_",".txt"),
		std::ios_base::out);
    njm::toFile(njm::toString(stats,"\n",""),
		njm::sett.datExt("sampStats_mean_"+file+"_",".txt"));
  }


  { // mode
    sObs.modelGen_r.mcmc.samples.setMode();
    std::vector<double> par = sObs.modelGen_r.mcmc.samples.getPar();
    sObs.modelGen_r.putPar(par.begin());

    /////////////
    // don't save mode to file (only save mean above)
    /////////////

    std::vector< std::vector<double> > stats;

    S s;
    Starts starts("startingLocations.txt");
    s.modelGen_r = sObs.modelGen_r;
    s.modelEst_r = s.modelGen_r;
    int r,t,R,T;
    R = numStats;
    T = sObs.sD.time;
    for(r = 0; r < R; ++r){
      s.modelGen_r.mcmc.samples.setRand();
      par = s.modelGen_r.mcmc.samples.getPar();
      s.modelGen_r.putPar(par.begin());
      s.modelEst_r.putPar(par.begin());

      s.reset(starts[r]);

      for(t = 0; t < T; ++t)
	s.nextPoint();
      h = s.sD.history;
      h.push_back(s.sD.status);

      stats.push_back(getStats(h,s.sD,s.fD));
    }

    njm::toFile(names,njm::sett.datExt("sampStats_mode_"+file+"_",".txt"),
		std::ios_base::out);
    njm::toFile(njm::toString(stats,"\n",""),
		njm::sett.datExt("sampStats_mode_"+file+"_",".txt"));
  }


  { // mle
    sObs.modelGen_r.setType(MLE);
    sObs.modelGen_r.fit(sObs.sD,sObs.tD,sObs.fD,sObs.dD,0);

    njm::toFile(njm::toString(sObs.modelGen_r.getPar()," ","\n"),
    	      njm::sett.datExt("sampSats_"+file+"_MLE_",".txt"));

    std::vector<double> par;

    std::vector< std::vector<double> > stats;

    S s;
    Starts starts("startingLocations.txt");
    s.modelGen_r = sObs.modelGen_r;
    s.modelEst_r = s.modelGen_r;
    int r,t,R,T;
    R = numStats;
    T = sObs.sD.time;
    for(r = 0; r < R; ++r){
      s.modelGen_r.mcmc.samples.setRand();
      par = s.modelGen_r.mcmc.samples.getPar();
      s.modelGen_r.putPar(par.begin());
      s.modelEst_r.putPar(par.begin());

      s.reset(starts[r]);

      for(t = 0; t < T; ++t)
	s.nextPoint();
      h = s.sD.history;
      h.push_back(s.sD.status);

      stats.push_back(getStats(h,s.sD,s.fD));
    }

    njm::toFile(names,njm::sett.datExt("sampStats_mle_"+file+"_",".txt"),
		std::ios_base::out);
    njm::toFile(njm::toString(stats,"\n",""),
		njm::sett.datExt("sampStats_mle_"+file+"_",".txt"));
  }



  // param estimates
  std::vector<std::vector<double> > parSamp;
  int i;
  for(i = 0; i < (numSamples-numBurn); ++i){
    sObs.modelGen_r.mcmc.samples.setPar(i);
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
			      " ","\n"),
  		njm::sett.datExt("sampStats_"+file+"_param_",".txt"),
  		std::ios_base::app);
  }

  for(i = 0; i < numBurn; ++i){
    sObs.modelGen_r.mcmc.samples.setPar(i,true);
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
			      " ","\n"),
  		njm::sett.datExt("sampStats_"+file+"_paramBurn_",".txt"),
  		std::ios_base::app);
  }

  // posterior mode
  sObs.modelGen_r.mcmc.samples.setMode();
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
  	      njm::sett.datExt("sampStats_"+file+"_paramMode_",".txt"));

  // posterior mean
  sObs.modelGen_r.mcmc.samples.setMean();
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
  	      njm::sett.datExt("sampStats_"+file+"_paramMean_",".txt"));

  // likelihood
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.ll,"\n",""),
  	      njm::sett.datExt("sampStats_"+file+"_ll_",".txt"));

  // likelihood at mean
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.llPt,"\n"),
  	      njm::sett.datExt("sampStats_"+file+"_llPt_",".txt"));

  // pD
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.pD,"\n"),
  	      njm::sett.datExt("sampStats_"+file+"_pD_",".txt"));

  // Dbar
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.Dbar,"\n"),
  	      njm::sett.datExt("sampStats_"+file+"_Dbar_",".txt"));

  // DIC
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.DIC,"\n"),
  	      njm::sett.datExt("sampStats_"+file+"_DIC_",".txt"));

}



int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  // int numSamples = 1000000, numBurn = 25000, numStats = 50000;
  // int numSamples = 50000, numBurn = 25000, numStats = 10000;
  int numSamples = 20000, numBurn = 10000, numStats = 10000;
  // int numSamples = 100, numBurn = 50, numStats = 50;
  // int numSamples = 10, numBurn = 5, numStats = 5;

#pragma omp parallel sections			\
  shared(numSamples,numBurn,numStats)
  {
#pragma omp section
    {
      runBayesP<ModelGravityGDist
		>("gravity",1,
		  numSamples,numBurn,numStats);
    }

#pragma omp section
    {
      runBayesP<Model2GravityGDist
		>("gravity2",0,
		  numSamples,numBurn,numStats);
    }

// #pragma omp section
//     {
//       runBayesP<ModelGravityGDistTrend
// 		>("gravityTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGravityGDistTrendPow
// 		>("gravityTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGravityGDistTrendPowCon
// 		>("gravityTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDist
// 		>("timeInf",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDistTrend
// 		>("timeInfTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDistTrendPow
// 		>("timeInfTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDistTrendPowCon
// 		>("timeInfTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDist
// 		>("timeInfExp",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDistTrend
// 		>("timeInfExpTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDistTrendPow
// 		>("timeInfExpTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDistTrendPowCon
// 		>("timeInfExpTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDist
// 		>("timeInfExpCaves",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDistTrend
// 		>("timeInfExpCavesTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDistTrendPow
// 		>("timeInfExpCavesTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDistTrendPowCon
// 		>("timeInfExpCavesTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelRad
// 		>("rad",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGDist
// 		>("gDist",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGDistPow
// 		>("gDistPow",0,
// 		  numSamples,numBurn,numStats);
//     }

  }

  // njm::sett.clean();
  return 0;
}
