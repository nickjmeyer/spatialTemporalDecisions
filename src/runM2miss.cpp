#include "runM2miss.hpp"

template <class ME, class PE>
void runM2(const std::string nameMod, const int i){
  typedef GravityTimeInfExpCavesModel MG;
  typedef GravityTimeInfExpCavesParam PG;
  
  typedef System<MG,PG,ME,PE> S;

  typedef ToyFeatures2<ME,PE> F;
  typedef FeaturesInt<F,ME,PE> FI;
  typedef RankAgent<F,ME,PE> RA;

  typedef M2QOptim<S,RA,FI,ME,PE> SPO;

  typedef OptimRunner<S,RA,SPO> R_RA;


  S s;

  s.modelGen.fitType = MLE;
  s.modelEst.fitType = MLE;

  RA ra;
  ra.name += "_" + nameMod + "_" + njm::toString(i,"",0,0);

  SPO spo;

  R_RA r_ra;
  

  int numReps = 100;

  njm::message("Priority Score: "
	       + njm::toString(r_ra.run(s,ra,spo,numReps,s.fD.finalT),""));
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MEexpcaves;
  typedef GravityTimeInfExpCavesParam PEexpcaves;
  
  // typedef RadiusModel MEradius;
  // typedef RadiusParam PEradius;

  runM2<MEexpcaves,PEexpcaves>("expcaves",2);
  // runM2<MEradius,PEradius>("radius",2);

  return 0;
}

