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

  M2NmOptim<System<GravityModel,GravityParam>,
  	    RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
  			 GravityModel,GravityParam>,
  	    FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
  			GravityModel,GravityParam>,
  	    GravityModel,GravityParam> m2;


  std::vector<double> sdStart,sdStop,sdDecay,sdJump;
  sdStart.push_back(1.0);
  sdStart.push_back(2.0);
  sdStart.push_back(5.0);

  sdStop.push_back(0.05);
  sdStop.push_back(0.1);
  sdStop.push_back(0.5);

  sdDecay.push_back(0.975);
  sdDecay.push_back(0.95);
  sdDecay.push_back(0.925);
  sdDecay.push_back(0.9);
    
  sdJump.push_back(5.0);
  sdJump.push_back(2.0);
  sdJump.push_back(1.5);
  sdJump.push_back(1.1);

  std::vector<double>::const_iterator start,stop,decay,jump;

  int i,n,r,N=12,R=50;

  std::vector<System<GravityModel,GravityParam> > sV;

  for(n=0; n<N; n++){

    s.reset();
    for(i=0; i<12; i++){
      if(i>=s.fD.trtStart)
	pA.applyTrt(s.sD,s.tD,s.fD,s.dD,s.model,s.estParam);

      s.updateStatus();
      s.nextPoint();
    }

    sV.push_back(s);
  }

  int threads0 = N, threads1 = omp_get_max_threads()/threads0;

  int m=0,M = N*R*sdStart.size()*sdStop.size()*sdDecay.size()*sdJump.size();

  int t,id,ncol=4+4+rA.f.numFeatures; // start, stop, decay, jump, id, n, q, t
  double q;
  int k;
  
  arma::mat results(0,ncol);

  omp_set_nested(1);
  
#pragma omp parallel for num_threads(threads0)			\
  shared(results,m,M,N,R,sdStart,sdStop,sdDecay,sdJump,ncol)	\
  firstprivate(sV,m2,rA)					\
  private(n,r,id,start,stop,decay,jump,t,q)
  for(n=0; n<N; n++){
    
#pragma omp parallel for num_threads(threads1)			\
  shared(results,n,m,M,N,R,sdStart,sdStop,sdDecay,sdJump,ncol)	\
  firstprivate(sV,m2,rA)					\
  private(r,id,start,stop,decay,jump,t,q)
    for(r=0; r<R; r++){
      id=0;
      for(start=sdStart.begin(); start!=sdStart.end(); start++){
	for(stop=sdStop.begin(); stop!=sdStop.end(); stop++){
	  for(decay=sdDecay.begin(); decay!=sdDecay.end(); decay++){
	    for(jump=sdJump.begin(); jump!=sdJump.end(); jump++){
	      rA.tp.weights.ones(rA.f.numFeatures);

	      m2.qEval.tp.sdStart = (*start);
	      m2.qEval.tp.sdStop = (*stop);
	      m2.qEval.tp.sdDecay = (*decay);
	      m2.qEval.tp.sdJump = (*jump);
	      
	      t = std::time(NULL);
	      m2.optim(sV.at(n),rA);
	      t = std::time(NULL) - t;

	      m2.qEval.preCompData(sV.at(n).sD,sV.at(n).fD);
	      m2.qEval.bellResFixData(sV.at(n).sD,sV.at(n).tD,sV.at(n).fD,
				      sV.at(n).dD,sV.at(n).model,
				      sV.at(n).estParam);
	      m2.qEval.bellResPolData(sV.at(n).sD.time,sV.at(n).fD,
				      sV.at(n).model,
				      sV.at(n).estParam,rA);
	      m2.qEval.solve();
	      q = m2.qEval.qFn(sV.at(n).sD,sV.at(n).tD,sV.at(n).fD,sV.at(n).dD,
			       sV.at(n).model,sV.at(n).estParam,rA);
	      
#pragma omp critical
	      {
		m++;
		results.resize(m,ncol);
		for(k=0; k<rA.f.numFeatures; k++)
		  results(m-1,k) = rA.tp.weights.at(k);
		
		results(m-1,k++) = (*start);
		results(m-1,k++) = (*stop);
		results(m-1,k++) = (*decay);
		results(m-1,k++) = (*jump);
		
		results(m-1,k++) = id++;
		results(m-1,k++) = n;
		results(m-1,k++) = q;
		results(m-1,k++) = t;
		
		results.save(njm::sett.datExt("results",".txt"),
			     arma::raw_ascii);

		printf("\rCompleted %08d out of %08d",m,M);
		fflush(stdout);
	      }
	    }
	  }
	}
      }
    }
  }
  printf("\n");


  results.save(njm::sett.datExt("results",".txt"),
	       arma::raw_ascii);
  

>>>>>>> m2OptimSA
  // njm::sett.clean();
  return 0;
}
