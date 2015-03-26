#include "runM1MLEupd.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef MG ME;
  typedef PG PE;

  typedef System<MG,PG,ME,PE> S;

  // typedef NoTrt<ME,PE> NT;
  // typedef ProximalAgent<ME,PE> PA;
  // typedef MyopicAgent<ME,PE> MA;
  
  typedef ToyFeatures2<ME,PE> F;
  typedef RankAgent<F,ME,PE> RA;

  typedef M1SpOptim<S,RA,ME,PE> SPO;

  // typedef VanillaRunner<S,NT> R_NT;
  // typedef VanillaRunner<S,PA> R_PA;
  // typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  // NT nt;
  // PA pa;
  // MA ma;
  RA ra;

  ra.tp.weights_r.zeros(ra.f.numFeatures);
  ra.tp.weights_r(2) = 1;
  
  SPO spo;
  // no tuning for right now....
  spo.tp.tune = 0;

  // R_NT r_nt;
  // R_PA r_pa;
  // R_MA r_ma;
  R_RA r_ra;
  

  int numReps = 100;
  

  // njm::message("  No treatment: "
  // 	       + njm::toString(r_nt.run(s,nt,numReps,s.fD.finalT),""));
  // njm::message("      Proximal: "
  // 	       + njm::toString(r_pa.run(s,pa,numReps,s.fD.finalT),""));
  // njm::message("        Myopic: "
  // 	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT),""));
  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT),""));

  return 0;
}

