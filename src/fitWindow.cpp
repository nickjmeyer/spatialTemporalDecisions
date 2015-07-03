#include "fitWindow.hpp"


template <class M>
void fitWindow(const std::string & ext,
	       const int numSamples, const int numBurn){
  njm::resetSeed(0);

  typedef System<M,M> S;

  S sObs("obsData_"+ext+".txt");

  sObs.modelGen_r.mcmc.load(sObs.sD.history,sObs.sD.status,sObs.fD);
  sObs.modelGen_r.mcmc.sample(numSamples,numBurn);

  // posterior mean
  sObs.modelGen_r.mcmc.samples.setMean();
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
	      njm::sett.datExt("postMean_"+ext+"_",".txt"));

  // mcmc samples
  int i;
  for(i = 0; i < (numSamples - numBurn); ++i){
    sObs.modelGen_r.mcmc.samples.setPar(i);
    if(i)
      njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
		  njm::sett.datExt("samples_"+ext+"_",".txt"),
		  std::ios_base::app);
    else
      njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
		  njm::sett.datExt("samples_"+ext+"_",".txt"),
		  std::ios_base::out);
  }

  // likelihood
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.ll,"\n",""),
  	      njm::sett.datExt("ll_"+ext+"_",".txt"));

  // likelihood at mean
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.llPt,"\n"),
  	      njm::sett.datExt("llPt_"+ext+"_",".txt"));

  // pD
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.pD,"\n"),
  	      njm::sett.datExt("pD_"+ext+"_",".txt"));

  // Dbar
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.Dbar,"\n"),
  	      njm::sett.datExt("Dbar_"+ext+"_",".txt"));

  // DIC
  njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.DIC,"\n"),
  	      njm::sett.datExt("DIC_"+ext+"_",".txt"));

  sObs.modelGen_r.setType(MLE);
  sObs.modelGen_r.fit(sObs.sD,sObs.tD,sObs.fD,sObs.dD,0);
  njm::toFile(njm::toString(sObs.modelGen_r.getPar()," ","\n"),
	      njm::sett.datExt("MLE_"+ext+"_",".txt"));
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelGravityGDist M;

  int i;
  // note can't put const variables in shared
  const int I = 7, win = 3, numSamples = 200, numBurn = 100;
#pragma omp parallel for num_threads(omp_get_num_threads())	\
  private(i)
  for(i = 0; i < I; ++i){
    std::string file = (njm::toString(i+1,"",0,0)
			+ "-"
			+ njm::toString(i+win,"",0,0));
    fitWindow<M>(file,numSamples,numBurn);
  }

  njm::sett.clean();
  return 0;
}
