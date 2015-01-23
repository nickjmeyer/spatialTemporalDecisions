#include "runSimEbola.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<EbolaModel,EbolaParam> s;
  
  NoTrt<EbolaModel,EbolaParam> nA;
  ProximalAgent<EbolaModel,EbolaParam> pA;
  MyopicAgent<EbolaModel,EbolaParam> mA;
  RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
	       EbolaModel,EbolaParam> rA;

  // no treatment
  VanillaRunner<System<EbolaModel,EbolaParam>,
		NoTrt<EbolaModel,EbolaParam> > nR;

  // proximal
  VanillaRunner<System<EbolaModel,EbolaParam>,
		ProximalAgent<EbolaModel,EbolaParam> > pR;
  
  // myopic
  FitOnlyRunner<System<EbolaModel,EbolaParam>,
		MyopicAgent<EbolaModel,EbolaParam> > mR;

  // rank agent fixed (1,1,1,1)
  FitOnlyRunner<System<EbolaModel,EbolaParam>,
		RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
			     EbolaModel,EbolaParam> > rR;

  // rank agent simple
  OptimRunner<System<EbolaModel,EbolaParam>,
	      RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
			   EbolaModel,EbolaParam>,
	      M1SimpleOptim<System<EbolaModel,EbolaParam>,
			    RankToyAgent<ToyFeatures1<EbolaModel,
						      EbolaParam>,
					 EbolaModel,EbolaParam>
			    > > rRSimple;
			    
  // rank agent sgd
  OptimRunner<System<EbolaModel,EbolaParam>,
	      RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
			   EbolaModel,EbolaParam>,
	      M1SgdOptim<System<EbolaModel,EbolaParam>,
			 RankToyAgent<ToyFeatures1<EbolaModel,
						   EbolaParam>,
				      EbolaModel,EbolaParam>
			 > > rRSgd;
    
  // rank agent hybrid
  OptimRunner<System<EbolaModel,EbolaParam>,
	      RankToyAgent<ToyFeatures1<EbolaModel,EbolaParam>,
			   EbolaModel,EbolaParam>,
	      M1HybridOptim<System<EbolaModel,EbolaParam>,
			    RankToyAgent<ToyFeatures1<EbolaModel,
						      EbolaParam>,
					 EbolaModel,EbolaParam>
			    > > rRHybrid;

  
  // m1 simple optim
  M1SimpleOptim<System<EbolaModel,EbolaParam>,
		RankToyAgent<ToyFeatures1<EbolaModel,
					  EbolaParam>,
			     EbolaModel,EbolaParam> > m1Simple;

  // m1 sgd optim
  M1SgdOptim<System<EbolaModel,EbolaParam>,
	     RankToyAgent<ToyFeatures1<EbolaModel,
				       EbolaParam>,
			  EbolaModel,EbolaParam> > m1Sgd;

  // m1 hybrid optim
  M1HybridOptim<System<EbolaModel,EbolaParam>,
		RankToyAgent<ToyFeatures1<EbolaModel,
					  EbolaParam>,
			     EbolaModel,EbolaParam> > m1Hybrid;


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
