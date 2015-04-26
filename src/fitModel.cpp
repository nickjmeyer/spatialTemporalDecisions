#include "fitModel.hpp"


template<class M>
void fitModel(const int numSamples, const int numBurn,
	      const int verbose, const std::string msg){
  typedef System<M,M> S;
  
  S s("obsData.txt");
  
  s.modelEst_r = s.modelGen_r;
  s.revert();

  s.modelGen.fitType = MLE;
  s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,0);

  std::vector<double> mlePar = s.modelGen.getPar();

  s.modelGen.mcmc.load(s.sD.history,s.sD.status,s.fD);
  s.modelGen.mcmc.sample(numSamples,numBurn);
  s.modelGen.mcmc.samples.setMean();
  std::vector<double> par = s.modelGen.mcmc.samples.getPar();
  s.modelGen.putPar(par.begin());

  // s.paramGen.save();

  std::vector<double> mcmcPar = s.modelGen.getPar();

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
    
      typedef ModelGravity M;

      fitModel<M>(numSamples,numBurn,1,msg);
    }

// #pragma omp section
//     {
//       std::string msg = "Gravity Model: linear time infected";
    
//       typedef ModelGravity M;

//       fitModel<M>(numSamples,numBurn,1,msg);
//     }

// #pragma omp section
//     {
//       std::string msg = "Gravity Model: exp time infected";
    
//       typedef ModelTime M;

//       fitModel<M>(numSamples,numBurn,1,msg);
//     }
    
// #pragma omp section
//     {
//       std::string msg = "Gravity Model: exp caves time infected";
    
//       typedef ModelTimeExpCaves M;

//       fitModel<M>(numSamples,numBurn,1,msg);
//     }
    
  }

  njm::sett.clean();
  return 0;
}
