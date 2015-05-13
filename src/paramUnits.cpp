#include "paramUnits.hpp"

void test(const std::string & name, const int cond){
  if(cond)
    printf("%32s: passed\n",name.c_str());
  else
    printf("%32s: failed ***\n",name.c_str());
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelGravityGDist MG;
  typedef MG ME;
  typedef System<MG,ME> S;

  S s;

  std::cout << "Testing ParamIntercept" << std::endl;
  ParamBase * pB = new ParamIntercept();
  pB->init(s.fD);

  std::vector<double> ans;
  ans = {0};
  test("init",pB->getPar() == ans);

  std::vector<double> pars = {1};
  std::vector<double>::const_iterator beg = pars.begin();

  beg = pB->putPar(beg);

  ans = {1};
  test("putPar",pB->getPar() == ans);

  test("putPar",beg == pars.end());


  std::vector<double> probs = {0,1,2};

  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);

  ans = {1,2,3};
  test("setFill",probs == ans);

  delete pB;

  std::cout << "Testing ParamBeta" << std::endl;
  
  pB = new ParamBeta();
  pB->init(s.fD);

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


  delete pB;


  pB = new ParamGravityGDist();
  pB->init(s.fD);
  std::cout << "Testing ParamGravityGDist" << std::endl;

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
      ans.push_back(-pars.at(0) * s.fD.gDist.at(i*s.fD.numNodes + j) /
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


  delete pB;




  pB = new ParamTrt();
  pB->init(s.fD);
  std::cout << "Testing ParamTrt" << std::endl;

  pars = {0,0};
  test("init",pB->getPar() == pars);

  pars = {-1,1};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  s.tD.a = s.tD.p = std::vector<int>(s.fD.numNodes,0);

  s.tD.a.at(0) = s.tD.a.at(3) = s.tD.a.at(7) = 1;
  s.tD.p.at(0) = s.tD.a.at(4) = s.tD.a.at(8) = 1;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(-pars.at(0) * double(s.tD.a.at(j)) -
		    pars.at(1) * double(s.tD.p.at(i)));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);
  
  s.tD.a.at(0) = s.tD.a.at(3) = 0;
  s.tD.p.at(0) = 0;

  s.tD.a.at(2) = s.tD.a.at(9) = 1;
  s.tD.p.at(5) = s.tD.p.at(6) = 1;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(-pars.at(0) * double(s.tD.a.at(j)) -
		    pars.at(1) * double(s.tD.p.at(i)));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);


  s.tD.a.at(0) = s.tD.a.at(3) = 1;

  s.tD.a.at(1) = s.tD.a.at(8) = 1;
  s.tD.p.at(2) = s.tD.p.at(7) = 1;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(-pars.at(0) * double(s.tD.a.at(j)) -
		    pars.at(1) * double(s.tD.p.at(i)));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill2",diff < 1e-10);

  delete pB;

  pB = new ParamTime();
  pB->init(s.fD);
  std::cout << "Testing ParamTime" << std::endl;

  pars = {0};
  test("init",pB->getPar() == pars);

  pars = {-1};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  s.sD.timeInf = std::vector<int>(s.fD.numNodes,0);
  s.sD.timeInf.at(0) = 2;
  s.sD.timeInf.at(1) = 7;
  s.sD.timeInf.at(4) = 8;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) * double(s.sD.timeInf.at(j)-1));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  s.sD.timeInf.at(1) = 4;
  s.sD.timeInf.at(3) = 3;
  s.sD.timeInf.at(4) = 10;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) * double(s.sD.timeInf.at(j)-1));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  s.sD.timeInf.at(1) = 1;
  s.sD.timeInf.at(7) = 3;
  s.sD.timeInf.at(4) = 0;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) * double(s.sD.timeInf.at(j)-1));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill2",diff < 1e-10);

  delete pB;
  

  pB = new ParamTimeExpCaves();
  pB->init(s.fD);
  std::cout << "Testing ParamTimeExpCaves" << std::endl;

  pars = {0};
  test("init",pB->getPar() == pars);

  pars = {-1};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  s.sD.timeInf = std::vector<int>(s.fD.numNodes,0);
  s.sD.timeInf.at(0) = 2;
  s.sD.timeInf.at(1) = 7;
  s.sD.timeInf.at(4) = 4;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) *
		    (std::exp(double(s.sD.timeInf.at(j)-1)/
			      s.fD.propCaves.at(j))-1.0));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  s.sD.timeInf.at(1) = 4;
  s.sD.timeInf.at(3) = 3;
  s.sD.timeInf.at(4) = 2;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) *
		    (std::exp(double(s.sD.timeInf.at(j)-1)/
			      s.fD.propCaves.at(j))-1.0));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  s.sD.timeInf.at(1) = 1;
  s.sD.timeInf.at(7) = 3;
  s.sD.timeInf.at(4) = 0;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0) *
		    (std::exp(double(s.sD.timeInf.at(j)-1)/
			      s.fD.propCaves.at(j)) - 1.0));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("modFill2",diff < 1e-10);

  
  delete pB;




  pB = new ParamRadius();
  pB->init(s.fD);
  std::cout << "Testing ParamRadius" << std::endl;

  pars = {0};
  test("init",pB->getPar() == pars);

  pars = {-1};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      if(s.fD.logGDist.at(i*s.fD.numNodes + j) > pars.at(0))
	ans.push_back(-500.0);
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  s.sD.timeInf.at(1) = 4;
  s.sD.timeInf.at(3) = 3;
  s.sD.timeInf.at(4) = 2;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      if(s.fD.logGDist.at(i*s.fD.numNodes + j) > pars.at(0))
	ans.push_back(-500.0);
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  s.sD.timeInf.at(1) = 1;
  s.sD.timeInf.at(7) = 3;
  s.sD.timeInf.at(4) = 0;

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      if(s.fD.logGDist.at(i*s.fD.numNodes + j) > pars.at(0))
	ans.push_back(-500.0);
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("modFill2",diff < 1e-10);

  
  delete pB;
  

  pB = new ParamGDistKern();
  pB->init(s.fD);
  std::cout << "Testing ParamGDistKern" << std::endl;

  pars = {0,0};
  test("init",pB->getPar() == pars);

  pars = {-1,1};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  double d;
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      d = s.fD.gDist.at(i*s.fD.numNodes + j);
      ans.push_back(-pars.at(0)*std::exp(-d*d/(2*std::exp(pars.at(1)))));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      d = s.fD.gDist.at(i*s.fD.numNodes + j);
      ans.push_back(-pars.at(0)*std::exp(-d*d/(2*std::exp(pars.at(1)))));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      d = s.fD.gDist.at(i*s.fD.numNodes + j);
      ans.push_back(-pars.at(0)*std::exp(-d*d/(2*std::exp(pars.at(1)))));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("modFill2",diff < 1e-10);

  
  delete pB;




  pB = new ParamTrend();
  pB->init(s.fD);
  std::cout << "Testing ParamTrend" << std::endl;

  pars = {0};
  test("init",pB->getPar() == pars);

  pars = {-1};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*double(s.sD.time));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*double(s.sD.time));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  s.sD.time += 3;
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*double(s.sD.time));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("modFill2",diff < 1e-10);

  
  delete pB;



  pB = new ParamTrendPow();
  pB->init(s.fD);
  std::cout << "Testing ParamTrendPow" << std::endl;

  pars = {0,0};
  test("init",pB->getPar() == pars);

  pars = {-1,-0.5};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*
		    std::pow(double(s.sD.time),pars.at(1)));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*
		    std::pow(double(s.sD.time),pars.at(1)));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  s.sD.time += 3;
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*
		    std::pow(double(s.sD.time),pars.at(1)));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("modFill2",diff < 1e-10);

  
  delete pB;

  

  
  pB = new ParamTrendPowCon();
  pB->init(s.fD);
  std::cout << "Testing ParamTrendPowCon" << std::endl;

  pars = {0,0};
  test("init",pB->getPar() == pars);

  pars = {-1,-0.5};
  beg = pars.begin();
  beg = pB->putPar(beg);
  test("putPar",pB->getPar() == pars);

  test("putPar return",beg == pars.end());

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*
		    std::pow(double(s.sD.time),-std::exp(pars.at(1))));
    }
  }

  probs = std::vector<double>(s.fD.numNodes*s.fD.numNodes,0.0);
  pB->setFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("setFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*
		    std::pow(double(s.sD.time),-std::exp(pars.at(1))));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }
    
  test("modFill",diff < 1e-10);

  ans.clear();
  ans.reserve(s.fD.numNodes*s.fD.numNodes);
  s.sD.time += 3;
  for(i = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j){
      ans.push_back(pars.at(0)*
		    std::pow(double(s.sD.time),-std::exp(pars.at(1))));
    }
  }

  pB->modFill(probs,s.sD,s.tD,s.fD,s.dD);
  
  diff = 0.0;
  for(i = 0, k = 0; i < s.fD.numNodes; ++i){
    for(j = 0; j < s.fD.numNodes; ++j, ++k){
      diff += (probs.at(k) - ans.at(k))*(probs.at(k) - ans.at(k));
    }
  }

  test("modFill2",diff < 1e-10);

  
  delete pB;
  

  
  njm::sett.clean();
  return 0;
}

