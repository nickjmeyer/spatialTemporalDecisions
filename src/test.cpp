#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDist MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  // typedef NoTrt<ME> NT;
  typedef ProxStocGDistAgent<ME> PA;
  // typedef MyopicAgent<ME> MA;
  
  // typedef ToyFeatures7<ME> F;
  // typedef RankAgent<F,ME> RA;
  // typedef OsspAgent<ME> OA;

  // typedef M1SpOptim<S,RA,ME> SPO;
  // typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  // typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunnerNS<S,PA> R_PA;
  // typedef VanillaRunnerNS<S,RA> R_RA;
  // typedef FitOnlyRunner<S,MA> R_MA;
  // typedef OptimRunner<S,RA,SPO> R_RA;
  // typedef OptimRunner<S,OA,OSSPO> R_OA;


  // S s("obsData.txt");
  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 50;
  Starts starts(numReps,s.fD.numNodes);
  // s.reset(starts[0]);

  // NT nt;
  PA pa;
  // MA ma;
  // RA ra;
  // ra.tp.jitterScale = std::numeric_limits<double>::max();
  // ra.tp.weights.ones();
  // ra.tp.weights *= -1.0;
  // OA oa;

  // SPO spo;
  // OSSPO osspo;

  // R_NT r_nt;
  R_PA r_pa;
  // R_MA r_ma;
  // R_RA r_ra;
  // R_OA r_oa;


  RunStats rs;

  int end = s.fD.finalT;

  omp_set_num_threads(1);

  rs = r_pa.run(s,pa,numReps,end,starts);
  njm::message("       Proximal: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  njm::sett.clean();
  return 0;
}
