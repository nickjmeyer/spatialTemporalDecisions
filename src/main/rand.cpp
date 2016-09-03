#include "rand.hpp"

RandParr randParr(omp_get_max_threads(),1000000);


RandParr::RandParr(const int numSource, const int numRand){
  this->numSource = numSource;
  this->numRand = numRand;
  
  
  vgRunif01.clear();
  vgRnorm01.clear();
  seeds.resize(numSource);
  
  int i;
  for(i = 0; i < numSource; ++i){
    vgRunif01.push_back(VGRunif01(boost::mt19937(i),
				  boost::uniform_real<>(0.0,1.0)));
    vgRnorm01.push_back(VGRnorm01(boost::mt19937(i),
				  boost::normal_distribution<>(0.0,1.0)));
    setSeed(i,i);
  }
  
  runif01Iter.resize(numSource);
  runif01End.resize(numSource);
  runif01Vals = std::vector<std::vector<double>
			    >(numSource,std::vector<double>(numRand,0));
  rnorm01Iter.resize(numSource);
  rnorm01End.resize(numSource);
  rnorm01Vals = std::vector<std::vector<double>
			    >(numSource,std::vector<double>(numRand,0));
  
  reset();
}


void RandParr::setSeed(const int source, const int seed){
  seeds.at(source) = seed;
}


void RandParr::reset(){
  int i;
  for(i = 0; i < numSource; ++i)
    reset(i);
}


void RandParr::reset(const int source){
  vgRunif01.at(source).engine().seed(seeds.at(source));
  vgRunif01.at(source).distribution().reset();
  
  vgRnorm01.at(source).engine().seed(seeds.at(source));
  vgRnorm01.at(source).distribution().reset();
  
  fillRunif01(source);
  fillRnorm01(source);
}


void RandParr::fillRunif01(const int source){
  std::vector<double>::iterator it,beg,end;
  
  beg = runif01Vals.at(source).begin();
  end = runif01Vals.at(source).end();

  runif01Iter.at(source) = beg;
  runif01End.at(source) = end;
  for(it = beg; it != end; ++it)
    *it = vgRunif01.at(source)(); // gen runif01 values
}


void RandParr::fillRnorm01(const int source){
  std::vector<double>::iterator it,beg,end;
  
  beg = rnorm01Vals.at(source).begin();
  end = rnorm01Vals.at(source).end();

  rnorm01Iter.at(source) = beg;
  rnorm01End.at(source) = end;
  for(it = beg; it != end; ++it)
    *it = vgRnorm01.at(source)(); // gen rnorm01 values
}


double RandParr::genRunif01(const int source){
  if(runif01Iter.at(source) == runif01End.at(source))
    fillRunif01(source);
  return *runif01Iter.at(source)++;
}


double RandParr::genRnorm01(const int source){
  if(rnorm01Iter.at(source) == rnorm01End.at(source))
    fillRnorm01(source);
  return *rnorm01Iter.at(source)++;
}


void njm::resetSeed(){
  randParr.reset(omp_get_thread_num());
}


void njm::resetSeed(const int seed){
  int source = omp_get_thread_num();
  randParr.setSeed(source,seed);
  randParr.reset(source);
}


void njm::resetSeedAll(){
  randParr.reset();
}


double njm::runif01(){
  return(randParr.genRunif01(omp_get_thread_num()));
}


double njm::runif(double a, double b){
  return njm::runif01()*(b-a)+a;
}


int njm::runifInterv(int min, int max){
  double r=njm::runif01();
  r*=(double)(max-min);
  r+=(double)min;
  return (int)floor(r);
}


int njm::rber(double p){
  double r=runif01();
  return (int)(r<p);
}


double njm::rnorm01(){
  return(randParr.genRnorm01(omp_get_thread_num()));
}


double njm::rgamma(double const alpha, double const beta){
  return beta*boost::math::gamma_p_inv(alpha,njm::runif01());
}


double njm::qnormal(double prob){
  return -std::sqrt(2.0)*boost::math::erfc_inv(2.0*prob);
}

double njm::pnormal(double quan){
  return .5*boost::math::erfc(-quan/std::sqrt(2.0));
}
