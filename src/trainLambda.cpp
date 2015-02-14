#include "trainLambda.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  int i,j,numSys=25;

  typedef GravityModel GM;
  typedef GravityParam GP;
  typedef System<GM,GP,GM,GP> S;
  typedef ToyFeatures2<GM,GP> F;
  typedef FeaturesInt<F,GM,GP> FI;
  typedef RankToyAgent<F,GM,GP> RA;
  typedef M2SaOptim<S,RA,FI,GM,GP> M2;
  
  std::vector<S> s(numSys);
  RA rA;  
  for(i=0; i<numSys; i++){
    s.at(i).paramEst_r = s.at(i).paramGen_r;
    s.at(i).reset();
    for(j=0; j<s.at(i).fD.finalT; j++){
      if(j>=s.at(i).fD.trtStart){
	rA.applyTrt(s.at(i).sD,s.at(i).tD,s.at(i).fD,s.at(i).dD,
		    s.at(i).modelEst,s.at(i).paramEst);
	s.at(i).updateStatus();
	s.at(i).nextPoint();
      }
    }
  }

  arma::mat results(0,3);
  
#pragma omp parallel num_threads(omp_get_max_threads())		\
  firstprivate(s,rA)						\
  shared(results,numSys)					\
  private(i,j)
  {
    M2 qO;

    int minLambda=0;
    double err,errOrr,minErr=-1;

    std::vector<double> lambdaVals;
    for(i=2; i<13; i++)
      for(j=1; j<10; j+=2)
	lambdaVals.push_back(j*std::pow(10,i));

	  
    int I = lambdaVals.size();

    int jp1;
    for(i=0; i<I; i++){
      qO.qEval.tp.lambda = lambdaVals.at(i);
      err=0;
      errOrr=0;
      for(j=0; j<numSys; j++){
	jp1 = (j+1 == numSys ? 0 : j+1);
	
	// train data set j
	qO.qEval.preCompData(s.at(j).sD,s.at(j).fD);
	qO.qEval.bellResFixData(s.at(j).sD,s.at(j).tD,s.at(j).fD,s.at(j).dD,
				s.at(j).modelEst,s.at(j).paramEst);
	qO.qEval.bellResPolData(s.at(j).sD.time,s.at(j).fD,
				s.at(j).modelEst,s.at(j).paramEst,rA);
	qO.qEval.solve();

	errOrr += qO.qEval.bellRes();

	// test data set j+1
	qO.qEval.preCompData(s.at(jp1).sD,s.at(jp1).fD);
	qO.qEval.bellResFixData(s.at(jp1).sD,s.at(jp1).tD,s.at(jp1).fD,
				s.at(jp1).dD,s.at(jp1).modelEst,
				s.at(jp1).paramEst);
	qO.qEval.bellResPolData(s.at(jp1).sD.time,s.at(jp1).fD,
				s.at(jp1).modelEst,s.at(jp1).paramEst,rA);

	// error
	err += qO.qEval.bellRes();
      }
      errOrr /= (double)numSys;
      err /= (double)numSys;
#pragma omp critical
      {
	results.resize(results.n_rows+1,results.n_cols);
	results(results.n_rows-1,0) = lambdaVals.at(i);
	results(results.n_rows-1,1) = errOrr;
	results(results.n_rows-1,2) = err;
	results.save(njm::sett.datExt("results",".txt"),arma::raw_ascii);
      }

      if(minErr < 0 || err < minErr){
	minErr = err;
	minLambda = i;
      }
    }


//     for(i=minLambda-3000; i<minLambda+3001; i+=100){
//       qO.qEval.tp.lambda = i;
//       err=0;
//       errOrr=0;
//       for(j=0; j<numSys; j++){
// 	jp1 = (j+1 == numSys ? 0 : j+1);

// 	// train data set j
// 	qO.qEval.preCompData(s.at(j).sD,s.at(j).fD);
// 	qO.qEval.bellResFixData(s.at(j).sD,s.at(j).tD,s.at(j).fD,s.at(j).dD,
// 				s.at(j).model,s.at(j).estParam);
// 	qO.qEval.bellResPolData(s.at(j).sD.time,s.at(j).fD,
// 				s.at(j).model,s.at(j).estParam,rA);
// 	qO.qEval.solve();
	
// 	errOrr += qO.qEval.bellRes();

// 	// test data set j+1
// 	qO.qEval.preCompData(s.at(jp1).sD,s.at(jp1).fD);
// 	qO.qEval.bellResFixData(s.at(jp1).sD,s.at(jp1).tD,s.at(jp1).fD,
// 				s.at(jp1).dD,s.at(jp1).model,
// 				s.at(jp1).estParam);
// 	qO.qEval.bellResPolData(s.at(jp1).sD.time,s.at(jp1).fD,
// 				s.at(jp1).model,s.at(jp1).estParam,rA);

// 	// error
// 	err += qO.qEval.bellRes();
//       }
//       err /= (double)numSys;
// #pragma omp critical
//       {
// 	results.resize(results.n_rows+1,results.n_cols);
// 	results(results.n_rows-1,0) = i;
// 	results(results.n_rows-1,1) = errOrr;
// 	results(results.n_rows-1,2) = err;
// 	results.save(njm::sett.datExt("results",".txt"),arma::raw_ascii);
//       }
//     }
  }  

  
  
  // njm::sett.clean();
  return 0;
}
