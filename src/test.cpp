#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDistTrendPowCon MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef ProximalGDistAgent<ME> PA;

  typedef NullOptim<S,PA,ME> NO;

  typedef IncremAgent<ME,PA,NO> IA;

  typedef FitOnlyRunner<S,IA> R_IA;
  

  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  IA ia;
  
  R_IA r_ia;

  RunStats rs;

  rs = r_ia.run(s,ia,numReps,s.fD.finalT,starts);
  njm::message("    IncremAgent: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  std::vector<double> par = s.modelGen_r.getPar();

  rs = r_nt.run(s,nt,numReps,s.fD.finalT,starts);
  njm::message(rs.smean());

  std::for_each(s.fD.gDist.begin(),s.fD.gDist.end(),
		[](double & x){
		  x += 1.0;
		});
  
  s.preCompData();
  
  s.modelGen_r = MG(s.fD);
  // s.modelGen_r.read();
  s.modelGen_r.putPar(par.begin());
  s.modelEst_r = ME(s.fD);
  s.modelGen_r.setType(INVALID);
  s.modelEst_r.setType(INVALID);

  rs = r_nt.run(s,nt,numReps,s.fD.finalT,starts);
  njm::message(rs.smean());

  std::for_each(s.fD.gDist.begin(),s.fD.gDist.end(),
		[](double & x){
		  x -= 1.0;
		});
  
  s.preCompData();
  
  s.modelGen_r = MG(s.fD);
  // s.modelGen_r.read();
  s.modelGen_r.putPar(par.begin());
  s.modelEst_r = ME(s.fD);
  s.modelGen_r.setType(INVALID);
  s.modelEst_r.setType(INVALID);

  rs = r_nt.run(s,nt,numReps,s.fD.finalT,starts);
  njm::message(rs.smean());

  njm::sett.clean();
  return 0;
}
