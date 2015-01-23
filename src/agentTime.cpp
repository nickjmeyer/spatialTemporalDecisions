#include "agentTime.hpp"
#include "omp.h"

using namespace std::chrono;

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam_r = s.genParam_r;
  RankAgent<GravityModel,GravityParam> rA;
  SweepAgent<GravityModel,GravityParam> sA;

  // N should match the number in the R code  
  int threads=omp_get_max_threads(),N=14,m=5,n0=2,R=2;
  threads=63;
  arma::mat timesRank(threads*R,N-m),timesSweep(threads*R,N-m);
  timesRank.zeros();
  timesSweep.zeros();


  int numNodesOrig=s.fD.numNodes;
  std::vector<double> distOrig;
  distOrig=s.fD.dist;
  std::vector<int> networkOrig;
  networkOrig=s.fD.network;

#pragma omp parallel num_threads(threads)			\
  shared(timesRank,timesSweep,N)				\
  firstprivate(s,rA,sA,distOrig,networkOrig,numNodesOrig)
  {
    double tick,tock;

    int i,j,k,numNodes,t=omp_get_thread_num();
    for(i=0; i<(N-m); i++){
      numNodes=(N-i)*(N-i)*n0*n0;
      
      // sD
      s.sD_r.time=0;
      
      s.sD_r.numInfected=1;
      s.sD_r.numNotInfec=numNodes-1;
      
      s.sD_r.notInfec.resize(s.sD_r.numNotInfec);

      s.sD_r.timeInf.resize(numNodes);
      s.sD_r.status.resize(numNodes);
      s.sD_r.history.clear();
      s.sD_r.history.push_back(s.sD_r.status);

      // tD
      s.tD_r.a.resize(numNodes);
      s.tD_r.p.resize(numNodes);

      s.tD_r.aPast.resize(numNodes);
      s.tD_r.pPast.resize(numNodes);

      // fD
      s.fD.numNodes=numNodes;

      s.fD.dist.clear();
      for(j=0; j<numNodes; j++){
	for(k=0; k<numNodes; k++){
	  s.fD.dist.push_back(distOrig.at(j*numNodesOrig + k));
	}
      }

      s.fD.caves.resize(numNodes);
      s.fD.covar.resize(numNodes*s.fD.numCovar);
      s.fD.fips.resize(numNodes);
      
      s.fD.network.clear();
      for(j=0; j<numNodes; j++){
	for(k=0; k<numNodes; k++){
	  s.fD.network.push_back(networkOrig.at(j*numNodesOrig + k));
	}
      }

      s.fD.centroidsLong.resize(numNodes);
      s.fD.centroidsLat.resize(numNodes);
      
      s.fD.subGraph.resize(numNodes);
      s.fD.betweenness.resize(numNodes);

      njm::message(i);

      for(j=0; j<R; j++){
	s.reset();
	while(s.value()<(1.0/3.0))
	  s.nextPoint();

	tick=duration_cast<nanoseconds>(high_resolution_clock::
					now().time_since_epoch()).count();
	rA.applyTrt(s.sD,s.tD,s.fD,s.dD,s.model,s.genParam);
	tock=duration_cast<nanoseconds>(high_resolution_clock::
					now().time_since_epoch()).count();
	std::fill(s.tD.a.begin(),s.tD.a.end(),0);
	std::fill(s.tD.p.begin(),s.tD.p.end(),0);

#pragma omp critical
	{
	  timesRank(t*R+j,(N-m)-i-1) = (tock-tick)*1.0e-9;
	  timesRank.save(njm::sett.datExt("timesRank",".txt"),
			 arma::raw_ascii);
	}

	tick=duration_cast<nanoseconds>(high_resolution_clock::
					now().time_since_epoch()).count();
	sA.applyTrt(s.sD,s.tD,s.fD,s.dD,s.model,s.genParam);
	tock=duration_cast<nanoseconds>(high_resolution_clock::
					now().time_since_epoch()).count();
	std::fill(s.tD.a.begin(),s.tD.a.end(),0);
	std::fill(s.tD.p.begin(),s.tD.p.end(),0);
#pragma omp critical
	{
	  timesSweep(t*R+j,(N-m)-i-1) = (tock-tick)*1.0e-9;
	  timesSweep.save(njm::sett.datExt("timesSweep",".txt"),
			  arma::raw_ascii);
	}
      }
      
    }
    
    
  }

  // njm::sett.clean();
  return 0;
}
