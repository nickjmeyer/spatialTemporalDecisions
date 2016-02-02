#include "testSgd.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<GravityModel,GravityParam> s;
  s.estParam_r = s.genParam_r;
  s.reset();
  
  RankToyAgent<GravityModel,GravityParam> rA;

  TestRunner<System,RankToyAgent,GravityModel,GravityParam,M1SgdOptim> rR;

  M1SgdOptim<System,RankToyAgent,GravityModel,GravityParam> m1Sgd;




  std::vector<double> aVals;
  // aVals.push_back(1);
  // aVals.push_back(3);
  // aVals.push_back(7);
  // aVals.push_back(10);
  // aVals.push_back(13);
  // aVals.push_back(17);
  // aVals.push_back(20);
  // aVals.push_back(40);
  aVals.push_back(10);
  aVals.push_back(20);
  aVals.push_back(40);
  aVals.push_back(50);
  aVals.push_back(75);
  aVals.push_back(100);
  aVals.push_back(150);
  int i;
  
  std::vector<std::pair<double,double> > abVals;
  for(i=0; i<(int)aVals.size(); i++)
    abVals.push_back(std::pair<double,double>(aVals.at(i),1));

  int numAbVals = abVals.size();
  
  int mcReps=300;
  
  double val=0;
  arma::mat results(0,5);
  for(i=0; i<numAbVals; i++){
    m1Sgd.tp.a=abVals.at(i).first;
    m1Sgd.tp.b=abVals.at(i).second;

//     // some momentum 50 reps
//     m1Sgd.tp.momRate=.5;
//     m1Sgd.tp.mcReps=50;

//     val = rR.run(s,rA,m1Sgd,mcReps,s.fD.finalT);

// #pragma omp critical
//     {
//       results.resize(results.n_rows+1,results.n_cols);
//       results(results.n_rows-1,0) = m1Sgd.tp.a;
//       results(results.n_rows-1,1) = m1Sgd.tp.b;
//       results(results.n_rows-1,2) = m1Sgd.tp.momRate;
//       results(results.n_rows-1,3) = m1Sgd.tp.mcReps;
//       results(results.n_rows-1,4) = val;
//       results.save(njm::sett.datExt("results",".txt"),arma::raw_ascii);
//     }

//     // no momentum 50 reps
//     m1Sgd.tp.momRate=0;
//     m1Sgd.tp.mcReps=50;

//     val = rR.run(s,rA,m1Sgd,mcReps,s.fD.finalT);

// #pragma omp critical
//     {
//       results.resize(results.n_rows+1,results.n_cols);
//       results(results.n_rows-1,0) = m1Sgd.tp.a;
//       results(results.n_rows-1,1) = m1Sgd.tp.b;
//       results(results.n_rows-1,2) = m1Sgd.tp.momRate;
//       results(results.n_rows-1,3) = m1Sgd.tp.mcReps;
//       results(results.n_rows-1,4) = val;
//       results.save(njm::sett.datExt("results",".txt"),arma::raw_ascii);
//     }


      
    // some momentum 100 reps
    m1Sgd.tp.momRate=.5;
    m1Sgd.tp.mcReps=100;

    val = rR.run(s,rA,m1Sgd,mcReps,s.fD.finalT);

#pragma omp critical
    {
      results.resize(results.n_rows+1,results.n_cols);
      results(results.n_rows-1,0) = m1Sgd.tp.a;
      results(results.n_rows-1,1) = m1Sgd.tp.b;
      results(results.n_rows-1,2) = m1Sgd.tp.momRate;
      results(results.n_rows-1,3) = m1Sgd.tp.mcReps;
      results(results.n_rows-1,4) = val;
      results.save(njm::sett.datExt("results",".txt"),arma::raw_ascii);
    }

//     // no momentum 100 reps
//     m1Sgd.tp.momRate=0;
//     m1Sgd.tp.mcReps=100;

//     val = rR.run(s,rA,m1Sgd,mcReps,s.fD.finalT);

// #pragma omp critical
//     {
//       results.resize(results.n_rows+1,results.n_cols);
//       results(results.n_rows-1,0) = m1Sgd.tp.a;
//       results(results.n_rows-1,1) = m1Sgd.tp.b;
//       results(results.n_rows-1,2) = m1Sgd.tp.momRate;
//       results(results.n_rows-1,3) = m1Sgd.tp.mcReps;
//       results(results.n_rows-1,4) = val;
//       results.save(njm::sett.datExt("results",".txt"),arma::raw_ascii);
//     }
  }
    
  
  // njm::sett.clean();
  
  return 0;
}
