#include "surface.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam = s.genParam;
  
  RankAgent<GravityModel,GravityParam> rA;

  PlainRunner<System,RankAgent,GravityModel,GravityParam> rR;


  std::vector< std::vector<double> > weightsAll;
  double r,x,y,z,w,theta,start=-.95,end=.95,by=.05,area,totArea,tol=.0000001;
  int Ntot=5000,Nmin=75,N,n;

  std::vector<double> weights;
  for(n=0; n<rA.numFeatures; n++){
    weights.resize(rA.numFeatures);
    std::fill(weights.begin(),weights.end(),0.0);
    
    weights.at(n) = 1;
    weightsAll.push_back(weights);

    weights.at(n) = -1;
    weightsAll.push_back(weights);
  }

  
  totArea=0;
  for(x=start; x<end+tol; x+=by)
    totArea += M_PI * (1-x*x);

  for(x=start; x<end+tol; x+=by){
    weights.push_back(x);
    area = M_PI*(1-x*x)/totArea;

    N = std::max(Nmin,((int)(((double)Ntot)*area))+1);
    for(n=0; n<N; n++){
      r = njm::runif01()*(1-x*x);
      theta = njm::runif01()*2*M_PI;

      y = std::sqrt(r)*std::cos(theta);
      z = std::sqrt(r)*std::sin(theta);

      w = std::sqrt(1.0 - x*x - y*y - z*z);
      
      weights.clear();
      weights.push_back(x);
      weights.push_back(y);
      weights.push_back(z);
      weights.push_back(w);
      weightsAll.push_back(weights);


      weights.clear();
      weights.push_back(x);
      weights.push_back(y);
      weights.push_back(z);
      weights.push_back(-w);
      weightsAll.push_back(weights);
    }
  }

  

  double value;

  int done=0,threads = (omp_get_max_threads() < 32 ? 1 : omp_get_max_threads());

  Ntot = weightsAll.size();
#pragma omp parallel for num_threads(threads)	\
  shared(Ntot,weightsAll,done)			\
  firstprivate(rR,s,rA)				\
  private(n,value)
  for(n=0; n<Ntot; n++){
#pragma omp critical
    {
      rA.tp.putPar(weightsAll.at(n));
    }
    
    value = rR.run(s,rA,150,s.fD.finalT);
    
#pragma omp critical
    {
      weightsAll.at(n).push_back(value);
      std::cout << "Completed " << std::setw(8) << ++done
		<< " out of " << std::setw(8) << Ntot
		<< "\r" << std::flush;
      if((done % 100) == 0)
	njm::toFile(njm::toString(weightsAll,"\n",""),
		    njm::sett.datExt("weights_",".txt"),
		    std::ios_base::out);
    }
  }
  njm::toFile(njm::toString(weightsAll,"\n",""),
	      njm::sett.datExt("weights_",".txt"),
	      std::ios_base::out);

  return 0;
}
