#include "runSimToy.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  
  NoTrt<GravityModel,GravityParam> nA;
  ProximalAgent<GravityModel,GravityParam> pA;
  MyopicAgent<GravityModel,GravityParam> mA;
  RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
	       GravityModel,GravityParam> rA;

  // no treatment
  VanillaRunner<System<GravityModel,GravityParam>,
		NoTrt<GravityModel,GravityParam> > nR;

  // proximal
  VanillaRunner<System<GravityModel,GravityParam>,
		ProximalAgent<GravityModel,GravityParam> > pR;
  
  // myopic
  FitOnlyRunner<System<GravityModel,GravityParam>,
		MyopicAgent<GravityModel,GravityParam> > mR;

  // rank agent fixed (1,1,1,1)
  FitOnlyRunner<System<GravityModel,GravityParam>,
		RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
			     GravityModel,GravityParam> > rR;

  // rank agent simple
  OptimRunner<System<GravityModel,GravityParam>,
	      RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
			   GravityModel,GravityParam>,
	      M1SimpleOptim<System<GravityModel,GravityParam>,
			    RankToyAgent<ToyFeatures2<GravityModel,
						      GravityParam>,
					 GravityModel,GravityParam>
			    > > rRSimple;
			    
  // rank agent sgd
  OptimRunner<System<GravityModel,GravityParam>,
	      RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
			   GravityModel,GravityParam>,
	      M1SgdOptim<System<GravityModel,GravityParam>,
			 RankToyAgent<ToyFeatures2<GravityModel,
						   GravityParam>,
				      GravityModel,GravityParam>
			 > > rRSgd;
    
  // rank agent hybrid
  OptimRunner<System<GravityModel,GravityParam>,
	      RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
			   GravityModel,GravityParam>,
	      M1HybridOptim<System<GravityModel,GravityParam>,
			    RankToyAgent<ToyFeatures2<GravityModel,
						      GravityParam>,
					 GravityModel,GravityParam>
			    > > rRHybrid;

  
  // m1 simple optim
  M1SimpleOptim<System<GravityModel,GravityParam>,
		RankToyAgent<ToyFeatures2<GravityModel,
					  GravityParam>,
			     GravityModel,GravityParam> > m1Simple;

  // m1 sgd optim
  M1SgdOptim<System<GravityModel,GravityParam>,
	     RankToyAgent<ToyFeatures2<GravityModel,
				       GravityParam>,
			  GravityModel,GravityParam> > m1Sgd;

  // m1 hybrid optim
  M1HybridOptim<System<GravityModel,GravityParam>,
		RankToyAgent<ToyFeatures2<GravityModel,
					  GravityParam>,
			     GravityModel,GravityParam> > m1Hybrid;


  int mcReps=300,numPoints = s.fD.finalT;
  njm::message("No Treatment");
  njm::message(nR.run(s,nA,mcReps,numPoints));
  
  njm::message("Proximal");
  njm::message(pR.run(s,pA,mcReps,numPoints));
  
  njm::message("Myopic");
  njm::message(mR.run(s,mA,mcReps,numPoints));
  
  njm::message("Rank");
  njm::message(rR.run(s,rA,mcReps,numPoints));
  
  njm::message("Rank Simple");
  njm::message(rRSimple.run(s,rA,m1Simple,mcReps,numPoints));

  njm::message("Rank Sgd");
  njm::message(rRSgd.run(s,rA,m1Sgd,mcReps,numPoints));
  
  njm::message("Rank Hybrid");
  njm::message(rRHybrid.run(s,rA,m1Hybrid,mcReps,numPoints));


  // njm::sett.clean();
  
  return 0;
}
