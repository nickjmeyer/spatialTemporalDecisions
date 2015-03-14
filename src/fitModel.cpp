#include "fitModel.hpp"


template<class M, class P>
void fitModel(const int numSamples, const int numBurn,
	      const int verbose, const std::string msg){
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

  // s.paramGen.save();

  std::vector<double> mcmcPar = s.paramGen.getPar();

#pragma omp critical
  {
    if(verbose > 0){
      std::cout << msg << std::endl;
      std::cout << " MLE: " + njm::toString(mlePar," ","\n");
      std::cout << "MCMC: " + njm::toString(mcmcPar," ","\n");
    }
  }
}

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  int numSamples=20,numBurn=10;

#pragma omp parallel sections			\
  shared(numSamples,numBurn)
  {
#pragma omp section
    {
      std::string msg = "Gravity Model: no time infected";
    
      typedef GravityModel M;
      typedef GravityParam P;

      fitModel<M,P>(numSamples,numBurn,1,msg);
    }

#pragma omp section
    {
      std::string msg = "Gravity Model: linear time infected";
    
      typedef GravityTimeInfModel M;
      typedef GravityTimeInfParam P;

      fitModel<M,P>(numSamples,numBurn,1,msg);
    }

#pragma omp section
    {
      std::string msg = "Gravity Model: exp time infected";
    
      typedef GravityTimeInfExpModel M;
      typedef GravityTimeInfExpParam P;

      fitModel<M,P>(numSamples,numBurn,1,msg);
    }
    
#pragma omp section
    {
      std::string msg = "Gravity Model: exp caves time infected";
    
      typedef GravityTimeInfExpCavesModel M;
      typedef GravityTimeInfExpCavesParam P;

      fitModel<M,P>(numSamples,numBurn,1,msg);
    }
    
  }

  njm::sett.clean();
  return 0;
}
