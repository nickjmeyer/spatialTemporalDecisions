#include "runM1.hpp"


template <class ME, class PE>
void runM1(const std::string nameMod, const int i){
  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef System<MG,PG,ME,PE> S;

  typedef ToyFeatures2<ME,PE> F;
  typedef RankAgent<F,ME,PE> RA;
  // typedef MyopicAgent<ME,PE> MA;
  
  typedef M1SpOptim<S,RA,ME,PE> SPO;

  // typedef FitOnlyRunner<S,MA> R_MA;  
  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;
  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  // MA ma;
  // ma.name += "_" + nameMod;
  RA ra;
  ra.name += "_" + nameMod + "_" + njm::toString(i,"",0,0);

  SPO spo;
  // no tuning for right now....
  spo.tp.tune = 0;

  // R_MA r_ma;  
  R_RA r_ra;
  

  int numReps = 100;
  
  
  // njm::message("        Myopic: "
  // 	       + njm::toString(r_ma.run(s,ma,numReps,s.fD.finalT),""));
  ra.tp.weights_r.zeros(ra.f.numFeatures);
  ra.tp.weights_r(i) = 1.0;
  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT),""));
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MEexpcaves;
  typedef GravityTimeInfExpCavesParam PEexpcaves;
  
  // typedef GravityTimeInfExpModel MEexp;
  // typedef GravityTimeInfExpParam PEexp;

  // typedef GravityTimeInfModel MElin;
  // typedef GravityTimeInfParam PElin;

  // typedef GravityModel MEgrav;
  // typedef GravityParam PEgrav;

  // typedef RangeModel MErange;
  // typedef RangeParam PErange;

  typedef RadiusModel MEradius;
  typedef RadiusParam PEradius;

  // typedef CaveModel MEcave;
  // typedef CaveParam PEcave;

  runM1<MEexpcaves,PEexpcaves>("expcaves",2);
  // runM1<MEexp,PEexp>("exp");
  // runM1<MElin,PElin>("lin");
  // runM1<MEgrav,PEgrav>("grav");
  // runM1<MErange,PErange>("range",i);
  runM1<MEradius,PEradius>("radius",2);
  // runM1<MEcave,PEcave>("caves",i);

  return 0;
}

