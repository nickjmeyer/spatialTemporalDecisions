#include "paramUnits.hpp"

void test(const std::string & name, const int cond){
  if(cond)
    printf("%32s: passed\n",name.c_str());
  else
    printf("%32s: failed ***\n",name.c_str());
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel MG;
  typedef MG ME;
  typedef System<MG,ME> S;

  S s;

  std::cout << "Testing ParamIntercept" << std::endl;
  ParamBase * pB = new ParamIntercept(s.fD);

  std::vector<double> ans;
  ans = {0};
  test("init",pB->getPar() == ans);

  std::vector<double> pars = {1};
  std::vector<double>::iterator beg = pars.begin();

  beg = pB->putPar(beg);

  ans = {1};
  test("putPar",pB->getPar() == ans);

  test("putPar",beg == pars.end());


  std::vector<double> probs = {0,1,2};

  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);

  ans = {1,2,3};
  test("setFill",probs == ans);

  pars = {2};
  pB->putPar(pars.begin());
  pB->updateFill(probs,s.sD,s.tD,s.fD,s.dD);

  ans = {2,3,4};
  test("updateFill",probs == ans);

  probs = {0,1,2};

  pars = {1};
  pB->putPar(pars.begin());
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);

  ans = {1,2,3};
  test("putPar,setFill",probs == ans);

  pars = {2};
  pB->putPar(pars.begin());
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  ans = {3,4,5};
  test("putPar,setFill",probs == ans);

  pB->updateFill(probs,s.sD,s.tD,s.fD,s.dD);

  ans = {4,5,6};
  test("updateFill",probs == ans);

  pars = {-1};
  pB->putPar(pars.begin());
  pB->updateFill(probs,s.sD,s.tD,s.fD,s.dD);

  ans = {1,2,3};
  test("putPar,updateFill",probs == ans);
  

  delete pB;

  std::cout << "Testing ParamBeta" << std::endl;
  
  pB = new ParamBeta(s.fD);

  ans = std::vector<double>(s.fD.numCovar,0);

  test("init",pB->getPar() == ans);

  pars = std::vector<double>(s.fD.numCovar,0);
  int i = 0;
  std::for_each(pars.begin(),pars.end(),
		[&i](double & x){ x = i++;});
  // pars is now {0,1,2,3,....numCovar - 1};

  beg = pars.begin();
  beg = pB->putPar(beg);

  ans = pars;
  test("putPar",pB->getPar() == ans);
  test("putPar return",beg == pars.end());

  
  int j;
  std::vector<double> covarBeta;
  covarBeta.resize(s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    covarBeta.at(i) = 0;
    for(j = 0; j < s.fD.numCovar; ++j){
      covarBeta.at(i) += s.fD.covar.at(i*s.fD.numCovar + j) * pars.at(j);
    }
  }
  ans = std::vector<double>(s.fD.numNodes * s.fD.numNodes,0);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.at(i*s.fD.numNodes + j) = covarBeta.at(i);
    }
  }

  probs = std::vector<double>(s.fD.numNodes * s.fD.numNodes,0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  double diff;
  diff = 0;
  for(i = 0; i < (s.fD.numNodes * s.fD.numNodes); ++i)
    diff += (ans.at(i) - probs.at(i))*(ans.at(i) - probs.at(i));

  test("setFill",diff < 1e-16);


  i = 0;
  std::for_each(pars.begin(),pars.end(),
		[&i](double & x){x = (i - 3)*(i + 2) - 1;});

  for(i = 0; i < s.fD.numNodes; ++i){
    covarBeta.at(i) = 0;
    for(j = 0; j < s.fD.numCovar; ++j){
      covarBeta.at(i) += s.fD.covar.at(i*s.fD.numCovar + j) * pars.at(j);
    }
  }
  ans = std::vector<double>(s.fD.numNodes * s.fD.numNodes,0);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.at(i*s.fD.numNodes + j) = covarBeta.at(i);
    }
  }

  pB->putPar(pars.begin());
  pB->updateFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0;
  for(i = 0; i < (s.fD.numNodes * s.fD.numNodes); ++i)
    diff += (ans.at(i) - probs.at(i))*(ans.at(i) - probs.at(i));

  test("updateFill",diff < 1e-16);

  delete pB;


  pB = new ParamGravity(s.fD);
  std::cout << "Testing ParamGravity" << std::endl;

  pars = {0,0};
  test("init",pB->getPar() == pars);

  pars = {1,1};
  beg = pars.begin();
  beg = pB->putPar(beg);

  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());


  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) * s.fD.dist.at(i*s.fD.numNodes + j) /
		    std::pow(s.fD.caves.at(i)*s.fD.caves.at(j),
			     pars.at(1)));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);

  diff = 0.0;
  int k;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);


  pars = {-2,3};
  beg = pars.begin();
  beg = pB->putPar(beg);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) * s.fD.dist.at(i*s.fD.numNodes + j) /
		    std::pow(s.fD.caves.at(i)*s.fD.caves.at(j),
			     pars.at(1)));
    }
  }

  pB->updateFill(probs,s.sD,s.tD,s.fD,s.dD);

  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("updateFill",diff < 1e-10);
  


  delete pB;
  
  
  njm::sett.clean();
  return 0;
}

