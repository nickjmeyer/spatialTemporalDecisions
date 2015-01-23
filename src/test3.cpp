#include "test3.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam_r = s.genParam_r;

  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
	       GravityModel,GravityParam> rA1;


  PlainRunner<System<GravityModel,GravityParam>,
	      RankToyAgent<ToyFeatures2<GravityModel,
					GravityParam>,
			   GravityModel,GravityParam> > tR1;

  resetRandomSeed();
  tR1.run(s,rA1,2000,s.fD.finalT);


  njm::message(s.fD.period);
  
  njm::sett.clean();
  return 0;
}
