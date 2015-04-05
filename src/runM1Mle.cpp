#include "runM1Mle.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef NoTrt<ME> NT;
  typedef ProximalAgent<ME> PA;
  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures2<ME> F;
  typedef RankAgent<F,ME> RA;

  typedef M1SpOptim<S,RA,ME> SPO;

  typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunner<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  NT nt;
  PA pa;
  MA ma;
  RA ra;

  ra.tp.weights_r.zeros(ra.f.numFeatures);
  ra.tp.weights_r(2) = 1;
  
  SPO spo;
  // no tuning for right now....
  spo.tp.tune = 0;

  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_RA r_ra;
  

  njm::message("  No treatment: "
	       + njm::toString(r_nt.run(s,nt,numReps,s.fD.finalT,starts),
			       ""));
  njm::message("      Proximal: "
	       + njm::toString(r_pa.run(s,pa,numReps,s.fD.finalT,starts),
			       ""));
  njm::message("        Myopic: "
	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT,starts),
			       ""));
  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts),
			       ""));

  return 0;
}

