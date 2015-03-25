#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

//   typedef GravityTimeInfExpCavesModel GM;
//   typedef GravityTimeInfExpCavesParam GP;
//   typedef GM EM;
//   typedef GP EP;

//   typedef System<GM,GP,EM,EP> S;

//   typedef ToyFeatures2<EM,EP> F;
//   typedef RankAgent<F,EM,EP> RA;

//   typedef FeaturesInt<F,EM,EP> FI;
//   typedef M2QOptim<S,RA,FI,EM,EP> OQ;

//   S s;

//   RA ra;
//   OQ oq;

//   s.modelGen.fitType = MLE;
//   s.modelEst.fitType = MLE;

//   s.paramEst_r = s.paramGen_r;
//   s.reset();

//   int i;
//   for(i = 0; i < 1; i++)
//     njm::runif01();
  
//   int t;
//   for(t = 0; t < 10; ++t){
//     if(t >= s.fD.trtStart){
//       ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.paramEst);
//     }
//     s.nextPoint();

//   }

//   std::cout << "value: " << s.value() << std::endl;
  
//   oq.qEval.preCompData(s.sD,s.fD);

//   oq.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst);

//   oq.qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,ra);

//   oq.qEval.buildRD();

// #pragma omp parallel for			\
//   private(i)					\
//   firstprivate(oq)
//   for(i = 0; i < 100; ++i){
//     std::cout << i << std::endl;
//     oq.qEval.solve();
//   }


  std::vector<int> ia,ja;
  std::vector<double> a,b,x;

  Eigen::SparseMatrix<double> mat(6,6);
  mat.insert(0,0) = 1;
  mat.insert(3,3) = 1;
  mat.insert(5,5) = 1;

  mat2Raw(mat,ia,ja,a);
  

  // njm::fromFile(ia,"ia_6.txt");
  // njm::fromFile(ja,"ja_6.txt");
  // njm::fromFile(a,"a_6.txt");
  // njm::fromFile(b,"b_6.txt");

  std::cout << "ia: [" << *std::min_element(ia.begin(),ia.end())
	    << ", " << *std::max_element(ia.begin(),ia.end()) << "]"
	    << "(" << ia.size() << ")" << std::endl;
  std::cout << "ja: [" << *std::min_element(ja.begin(),ja.end())
	    << ", " << *std::max_element(ja.begin(),ja.end()) << "]"
	    << "(" << ja.size() << ")" << std::endl;
  std::cout << "a: [" << *std::min_element(a.begin(),a.end())
	    << ", " << *std::max_element(a.begin(),a.end()) << "]"
	    << "(" << a.size() << ")" << std::endl;
  std::cout << "b: [" << *std::min_element(b.begin(),b.end())
	    << ", " << *std::max_element(b.begin(),b.end()) << "]"
	    << "(" << b.size() << ")" << std::endl;


  x = pardisoSymWrap(ia,ja,a,b);

  // ia = {1,2,3};
  // ja = {1,2};
  // a = {1,0};
  // b = {1,1};

  // x = pardisoSymWrap(ia,ja,a,b);  


  
  njm::sett.clean();
  return 0;
}
