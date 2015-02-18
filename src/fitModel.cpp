#include "fitModel.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef GravityModel EM;
  typedef GravityParam EP;

  typedef System<GM,GP,EM,EP> S;
  
  S s("dataObs.txt");
  
  s.paramEst_r = s.paramGen_r;
  s.reset();

  s.modelGen.fit(s.sD,s.tD,s.fD,s.dD,s.paramGen);

  s.paramGen.save();

  njm::message(s.paramGen.getPar());


  njm::sett.clean();
  return 0;
}
