#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam,
	 GravityModel,GravityParam> s;
  s.paramEst_r = s.paramGen_r;
  s.reset();

  ProximalAgent<GravityModel,GravityParam> pA;
  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
	       GravityModel,GravityParam> rA;

  M2NmOptim<System<GravityModel,GravityParam,
		   GravityModel,GravityParam>,
  	    RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
  			 GravityModel,GravityParam>,
  	    FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
  			GravityModel,GravityParam>,
  	    GravityModel,GravityParam> m2;


  njm::sett.clean();
  return 0;
}
