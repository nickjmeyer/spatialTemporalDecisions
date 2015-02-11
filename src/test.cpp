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

<<<<<<< HEAD
  M2NmOptim<System<GravityModel,GravityParam,
		   GravityModel,GravityParam>,
=======
  M2NmOptim<System<GravityModel,GravityParam>,
>>>>>>> m2OptimRand
  	    RankToyAgent<ToyFeatures2<GravityModel,GravityParam>,
  			 GravityModel,GravityParam>,
  	    FeaturesInt<ToyFeatures2<GravityModel,GravityParam>,
  			GravityModel,GravityParam>,
  	    GravityModel,GravityParam> m2;

<<<<<<< HEAD
=======
  std::vector<int> sampVals;
  sampVals.push_back(10);
  sampVals.push_back(50);
  sampVals.push_back(100);
  sampVals.push_back(250);
  sampVals.push_back(500);
  sampVals.push_back(1000);
  sampVals.push_back(10000);

  std::vector<int>::const_iterator samp;

  int i,n,r,N=12,R=20;

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

  int threads0 = N, threads1 = 4;
>>>>>>> m2OptimRand

  int m=0,M = N*R*sampVals.size();

  int t,id,ncol=1+4+rA.f.numFeatures; // start, stop, decay, jump, id, n, q, t
  double q;
  int k;
  
  arma::mat results(0,ncol);

  omp_set_nested(1);
  
#pragma omp parallel for num_threads(threads0)			\
  shared(results,m,M,N,R,sampVals,ncol)			\
  firstprivate(sV,m2,rA)					\
  private(n,r,id,samp,t,q)
  for(n=0; n<N; n++){
    
#pragma omp parallel for num_threads(threads1)			\
  shared(results,n,m,M,N,R,sampVals,ncol)	\
  firstprivate(sV,m2,rA)					\
  private(r,id,samp,t,q)
    for(r=0; r<R; r++){
      id=0;
      for(samp=sampVals.begin(); samp!=sampVals.end(); samp++){
	rA.tp.weights.ones(rA.f.numFeatures);

	m2.qEval.tp.numSamp = (*samp);
	      
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
		
	  results(m-1,k++) = (*samp);
		
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
  printf("\n");


  results.save(njm::sett.datExt("results",".txt"),
	       arma::raw_ascii);

  
  // njm::sett.clean();
  return 0;
}
