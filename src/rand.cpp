#include "rand.hpp"

// static unsigned randomSeed=std::time(NULL);
#ifdef RANDOM_SEED__
static unsigned randomSeed= RANDOM_SEED__;
#else
static unsigned randomSeed= 8;
#endif
static boost::mt19937 rng( randomSeed );
static boost::uniform_real<> unif01(0.0,1.0);
static boost::normal_distribution<> norm01(0.0, 1.0);
static boost::variate_generator<boost::mt19937,
				boost::uniform_real<>
				> genUnif01(rng, unif01);
static boost::variate_generator<boost::mt19937,
				boost::normal_distribution<> 
				> genNorm01(rng, norm01);

static RandParr randParr(1000000);



RandParr::RandParr(int numRand_){
  numRand=numRand_;
  initialize();
}


void RandParr::reset(){
  genUnif01.engine().seed(randomSeed);
  genUnif01.distribution().reset();
  
  genNorm01.engine().seed(randomSeed);
  genNorm01.distribution().reset();
  
  initialize();
}


void resetRandomSeed(){
  randParr.reset();
}



void RandParr::initialize(){
  int numThreads=omp_get_max_threads();
  // int numThreads=64;
  int i,j;
  numRunif01=numRand;
  numRnorm01=numRand;
  std::vector<double> threadRunif01(numRand), threadRunif01_fixed(numRand);
  std::vector<double> threadRnorm01(numRand), threadRnorm01_fixed(numRand);

  // allocate containers sizes
  runif01Vals.resize(numThreads);
  runif01Vals_fixed.resize(numThreads);
  rnorm01Vals.resize(numThreads);
  rnorm01Vals_fixed.resize(numThreads);
  
  // sample values for each thread
  for(i=0; i<numThreads; i++){
    for(j=0; j<numRand; j++){
      threadRunif01.at(j)=genUnif01();
      threadRunif01_fixed.at(j)=genUnif01();
      threadRnorm01.at(j)=genNorm01();
      threadRnorm01_fixed.at(j)=genNorm01();
    }
    // load the samples into the thread slots
    runif01Vals.at(i) = threadRunif01;
    runif01Vals_fixed.at(i) = threadRunif01_fixed;

    rnorm01Vals.at(i) = threadRnorm01;
    rnorm01Vals_fixed.at(i) = threadRnorm01_fixed;
  }

  // store iterators
  runif01Iter.clear();
  runif01End.clear();
  rnorm01Iter.clear();
  rnorm01End.clear();
  isfixed.clear();
  for(i=0; i<numThreads; i++){
    runif01Iter.push_back(runif01Vals.at(i).begin());
    runif01End.push_back(runif01Vals.at(i).end());
    
    rnorm01Iter.push_back(rnorm01Vals.at(i).begin());
    rnorm01End.push_back(rnorm01Vals.at(i).end());

    isfixed.push_back(false);
  }

  runif01Iter_hold.resize(numThreads);
  runif01End_hold.resize(numThreads);
  rnorm01Iter_hold.resize(numThreads);
  rnorm01End_hold.resize(numThreads);
}


void RandParr::fixSeed(const int fix){
  int thread=omp_get_thread_num();
  isfixed.at(thread)=fix;
  if(fix){
    runif01Iter_hold.at(thread) = runif01Iter.at(thread);
    rnorm01Iter_hold.at(thread) = rnorm01Iter.at(thread);
    
    runif01Iter.at(thread) = runif01Vals_fixed.at(thread).begin();
    rnorm01Iter.at(thread) = rnorm01Vals_fixed.at(thread).begin();
  }
  else{
    runif01Iter.at(thread) = runif01Iter_hold.at(thread);
    rnorm01Iter.at(thread) = rnorm01Iter_hold.at(thread);
  }
}


void fixRandomSeed(const int fix){
  randParr.fixSeed(fix);
}



double RandParr::getRunif01(){
  int thread=omp_get_thread_num();
  double r=*(runif01Iter.at(thread)++);

  // if there are no more samples, refill
  if(runif01Iter.at(thread) == runif01End.at(thread) && !isfixed.at(thread)){
#pragma omp critical
    {
      int i;
      for(i=0; i<numRunif01; i++)
	runif01Vals.at(thread).at(i)=genUnif01();
      runif01Iter.at(thread)=runif01Vals.at(thread).begin();
    }
  }
  if(runif01Iter.at(thread) == runif01End.at(thread) && isfixed.at(thread)){
#pragma omp critical
    {
      runif01Iter.at(thread) = runif01Vals_fixed.at(thread).begin();
    }
  }
  return r;
}


double RandParr::getRnorm01(){
  int thread=omp_get_thread_num();
  double r=*(rnorm01Iter.at(thread)++);

  // if there are no more samples, refill
  if(rnorm01Iter.at(thread) == rnorm01End.at(thread) && !isfixed.at(thread)){
#pragma omp critical
    {
      int i;
      for(i=0; i<numRnorm01; i++)
	rnorm01Vals.at(thread).at(i)=genNorm01();
      rnorm01Iter.at(thread)=rnorm01Vals.at(thread).begin();
    }
  }
  if(rnorm01Iter.at(thread) == rnorm01End.at(thread) && isfixed.at(thread)){
#pragma omp critical
    {
      rnorm01Iter.at(thread) = rnorm01Vals_fixed.at(thread).begin();
    }
  }
  return r;
}


double njm::runif01(){
  return(randParr.getRunif01());
}


double njm::runif(double a, double b){
  return runif01()*(b-a)+a;
}


int njm::runifInterv(int min, int max){
  double r=runif01();
  r*=(double)(max-min);
  r+=(double)min;
  return (int)floor(r);
}


int njm::rber(double p){
  double r=runif01();
  return (int)(r<p);
}


double njm::rnorm01(){
  return(randParr.getRnorm01());
}


double njm::rgamma(double const alpha, double const beta){
  return beta*boost::math::gamma_p_inv(alpha,runif01());
}


double njm::qnormal(double prob){
  return -std::sqrt(2.0)*boost::math::erfc_inv(2.0*prob);
}

double njm::pnormal(double quan){
  return .5*boost::math::erfc(-quan/std::sqrt(2.0));
}
