#include "fitModel.hpp"


template<class M, class P>
void fitModel(const int numSamples,const int numBurn,const int verbose){
  typedef System<M,P,M,P> S;
  
  S s("obsData.txt");
  
  s.paramEst_r = s.paramGen_r;
  s.reset();

  s.modelGen.fitType = MLE;
  s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramGen);

  std::vector<double> mlePar = s.paramGen.getPar();

  s.modelGen.mcmc.load(s.sD.history,s.sD.status,s.fD);
  s.modelGen.mcmc.sample(20000,10000);
  s.modelGen.mcmc.samples.setMean();
  s.paramGen.putPar(s.modelGen.mcmc.samples.getPar());

  s.paramGen.save();

  std::vector<double> mcmcPar = s.paramGen.getPar();

#pragma omp critical (verbose)
  {
    if(verbose){
      njm::message("MCMC: " + njm::toString(mlePar," "));
      njm::message("MCMC: " + njm::toString(mcmcPar," "));
    }
  }
}

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  int numSamples=20000,numBurn=5000;

#pragma omp parallel sections			\
  shared(numSamples,numBurn)
  {
#pragma omp section
    {
      njm::message("Gravity Model: no time infected");
    
      typedef GravityModel M;
      typedef GravityParam P;

      fitModel<M,P>(numSamples,numBurn,1);
    }
    
  }

  njm::sett.clean();
  return 0;
}
