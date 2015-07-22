#include "fitWindow.hpp"


template <class M>
void fitWindow(const std::string & ext,
	       const int numSamples, const int numBurn){
  njm::resetSeed(0);

  typedef System<M,M> S;

  S sObs("obsData_"+ext+".txt");

  std::string name = sObs.modelGen_r.name;

  sObs.modelGen_r.mcmc.load(sObs.sD.history,sObs.sD.status,sObs.fD);
  sObs.modelGen_r.mcmc.sample(numSamples,numBurn);

  // posterior mean
  sObs.modelGen_r.mcmc.samples.setMean();
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
  	      njm::sett.datExt(name + "_postMean_"+ext+"_",".txt"));

  std::vector<double> par = sObs.modelGen_r.mcmc.samples.getPar();

  // posterior mode
  sObs.modelGen_r.mcmc.samples.setMode();
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
  	      njm::sett.datExt(name + "_postMode_"+ext+"_",".txt"));



  // mcmc samples
  int i;
  for(i = 0; i < (numSamples - numBurn); ++i){
    sObs.modelGen_r.mcmc.samples.setPar(i);
    if(i)
      njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
  				" ","\n"),
  		  njm::sett.datExt(name+"_samples_"+ext+"_",".txt"),
  		  std::ios_base::app);
    else
      njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
  				" ","\n"),
  		  njm::sett.datExt(name+"_samples_"+ext+"_",".txt"),
  		  std::ios_base::out);
  }

  // likelihood
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.ll,"\n",""),
  	      njm::sett.datExt(name+"_ll_"+ext+"_",".txt"));

  // likelihood at mean
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.llPt,"\n"),
  	      njm::sett.datExt(name+"_llPt_"+ext+"_",".txt"));

  // pD
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.pD,"\n"),
  	      njm::sett.datExt(name+"_pD_"+ext+"_",".txt"));

  // Dbar
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.Dbar,"\n"),
  	      njm::sett.datExt(name+"_Dbar_"+ext+"_",".txt"));

  // DIC
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.DIC,"\n"),
  	      njm::sett.datExt(name+"_DIC_"+ext+"_",".txt"));

  sObs.modelGen_r.setType(MLE);
  sObs.modelGen_r.fit(sObs.sD,sObs.tD,sObs.fD,sObs.dD,0);
  njm::toFile(njm::toString(sObs.modelGen_r.getPar()," ","\n"),
	      njm::sett.datExt(name+"_MLE_"+ext+"_",".txt"));





  typedef NoTrt<M> NT;
  typedef VanillaRunnerNS<S,NT> R;

  S s;
  s.modelGen_r.putPar(par.begin());
  s.modelEst_r = s.modelGen_r;

  NT nt;

  R r;
  Starts starts("startingLocations.txt");
  RunStats rs;

  rs = r.run(s,nt,500,s.fD.finalT,starts);
  njm::toFile(njm::toString(rs.smean(),"\n"),
	      njm::sett.datExt(name+"_fcReg_"+ext+"_",".txt"));

  s.fD.forecastFlat = true;

  rs = r.run(s,nt,500,s.fD.finalT,starts);
  njm::toFile(njm::toString(rs.smean(),"\n"),
	      njm::sett.datExt(name+"_fcFlat_"+ext+"_",".txt"));
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  omp_set_nested(1);

#pragma omp parallel sections
  {
#pragma omp section
    {
      typedef ModelGravityGDist M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }

#pragma omp section
    {
      typedef ModelTimeExpCavesGDistTrendPowCon M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }

#pragma omp section
    {
      typedef ModelTimeExpGDistTrendPow M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }


#pragma omp section
    {
      typedef ModelTimeGDistTrendPow M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }

#pragma omp section
    {
      typedef ModelTimeExpCavesGDistTrend M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }

#pragma omp section
    {
      typedef ModelTimeExpGDistTrend M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }

#pragma omp section
    {
      typedef ModelTimeGDistTrend M;

      int i;
      // note can't put const variables in shared
      const int I = 6, win = 3, numSamples = 20000, numBurn = 10000;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
  private(i)
      for(i = 0; i < I; ++i){
	std::string file = (njm::toString(i,"",0,0)
			    + "-"
			    + njm::toString(i+win-1,"",0,0));
	fitWindow<M>(file,numSamples,numBurn);
      }
    }
  }
  return 0;
}
