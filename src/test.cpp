#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDist MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  // typedef NoTrt<ME> NT;
  // typedef ProximalGDistAgent<ME> PA;
  // typedef MyopicAgent<ME> MA;
  
  // typedef ToyFeatures4<ME> F;
  // typedef RankAgent<F,ME> RA;
  // typedef OsspAgent<ME> OA;

  // typedef M1SpOptim<S,RA,ME> SPO;
  // typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  // typedef VanillaRunner<S,NT> R_NT;
  // typedef VanillaRunner<S,PA> R_PA;
  // typedef FitOnlyRunner<S,MA> R_MA;
  // typedef OptimRunner<S,RA,SPO> R_RA;
  // typedef OptimRunner<S,OA,OSSPO> R_OA;


  S s("obsData.txt");
  // S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  // int numReps = 96;
  Starts starts("startingLocations.txt");
  s.reset(starts[0]);

  // NT nt;
  // PA pa;
  // MA ma;
  // RA ra;
  // OA oa;

  // SPO spo;
  // OSSPO osspo;

  // R_NT r_nt;
  // R_PA r_pa;
  // R_MA r_ma;
  // R_RA r_ra;
  // R_OA r_oa;


  // RunStats rs;

  // rs = r_nt.run(s,nt,numReps,s.fD.finalT,starts);
  // njm::message("   No treatment: "
  // 	       + njm::toString(rs.smean(),"")
  // 	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // rs = r_pa.run(s,pa,numReps,s.fD.finalT,starts);
  // njm::message("       Proximal: "
  // 	       + njm::toString(rs.smean(),"")
  // 	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // rs = r_ma.run(s,ma,numReps,s.fD.finalT,starts);
  // njm::message("         Myopic: "
  // 	       + njm::toString(rs.smean(),"")
  // 	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // rs = r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts);
  // njm::message("  Policy Search: "
  // 	       + njm::toString(rs.smean(),"")
  // 	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
  // njm::message("One Step Polish: "
  // 	       + njm::toString(rs.smean(),"")
  // 	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  // int i;
  // for(i = 0; i < 8; ++i){
  //   s.updateStatus();
  //   s.nextPoint();
  // }

  std::cout << "value: " << njm::toString(s.value(),"\n");

  int numSamples = 500, numBurn = 100;

#pragma omp parallel sections
  {
#pragma omp section
    {
      GravityMcmc mcmc;
      mcmc.load(s.sD.history,s.sD.status,s.fD);
      mcmc.sample(numSamples,numBurn);

      mcmc.samples.setMean();
      std::cout << "Gravity::par: " << std::endl
		<< njm::toString(mcmc.samples.getPar()," ","")
		<< std::endl;;
    }
    
#pragma omp section
    {
      GravityTimeInfMcmc mcmc;
      mcmc.load(s.sD.history,s.sD.status,s.fD);
      mcmc.sample(numSamples,numBurn);

      mcmc.samples.setMean();
      std::cout << "GravityTimeInf::par: " << std::endl
		<< njm::toString(mcmc.samples.getPar()," ","")
		<< std::endl;;
    }
    
#pragma omp section
    {
      GravityTimeInfExpMcmc mcmc;
      mcmc.load(s.sD.history,s.sD.status,s.fD);
      mcmc.sample(numSamples,numBurn);

      mcmc.samples.setMean();
      std::cout << "GravityTimeInfExp::par: " << std::endl
		<< njm::toString(mcmc.samples.getPar()," ","")
		<< std::endl;
    }  

#pragma omp section
    {
      GravityTimeInfExpCavesMcmc mcmc;
      mcmc.load(s.sD.history,s.sD.status,s.fD);
      mcmc.sample(numSamples,numBurn);

      mcmc.samples.setMean();
      std::cout << "GravityTimeInfExpCaves::par: " << std::endl
		<< njm::toString(mcmc.samples.getPar()," ","")
		<< std::endl;;
    }
  }

  njm::sett.clean();
  return 0;
}
