#include "test.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfExpCavesModel GM;
  typedef GravityTimeInfExpCavesParam GP;
  typedef GM EM;
  typedef GP EP;

  typedef System<GM,GP,EM,EP> S;

  typedef ToyFeatures2<EM,EP> F;
  typedef RankAgent<F,EM,EP> RA;

  typedef FeaturesInt<F,EM,EP> FI;
  typedef M2QOptim<S,RA,FI,EM,EP> OQ;

  S s;

  s.modelGen.infProbs(s.sD,s.tD,s.fD,s.dD,s.paramGen);

  double prob = 0;
  for(int i = 0; i < s.sD.numNotInfec; ++i)
    prob += s.paramGen.infProbs.at(i);

  std::cout << prob/double(s.sD.numNotInfec) << std::endl;

  
  // RA ra;
  // OQ oq;

  // s.modelGen.fitType = MLE;
  // s.modelEst.fitType = MLE;

  // s.paramEst_r = s.paramGen_r;
  // s.reset();

  // int i;
  // for(i = 0; i < 0; i++)
  //   njm::runif01();
  
  // int t;
  // for(t = 0; t < s.fD.finalT; ++t){
  //   if(t >= s.fD.trtStart){
  //     ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.paramEst);
  //   }
  //   s.nextPoint();

  // }

  // std::cout << "value: " << s.value() << std::endl;
  
  // oq.qEval.preCompData(s.sD,s.fD);

  // oq.qEval.bellResFixData(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst);

  // oq.qEval.bellResPolData(s.sD.time,s.fD,s.modelEst,s.paramEst,ra);

  // oq.qEval.buildRD();

  // std::cout << oq.qEval.R.sum() << " >> "
  // 	    << oq.qEval.D0.sum() << " >> "
  // 	    << oq.qEval.D1.sum() << " >> "
  // 	    << oq.qEval.D.sum()
  // 	    << std::endl;

  // oq.qEval.solve();

  // std::cout << oq.qEval.qFn(s.sD,s.tD,s.fD,s.dD,s.modelEst,s.paramEst,ra)
  // 	    << " >>>>> " << oq.qEval.bellRes()
  // 	    << std::endl;
  
  // std::cout << "lambda before: " << oq.qEval.tp.lambda << std::endl;
  // oq.qEval.tune(s.sD.status);
  // std::cout << " lambda after: " << oq.qEval.tp.lambda << std::endl;

  // std::cout << "optimizing...."
  // 	    << std::endl;

  // oq.qEval.tp.lambda = 2.76e+06;
  // oq.optim(s,ra);


  // std::vector<int> nodes;
  // for(i = 0; i < s.fD.numNodes; ++i)
  //   nodes.push_back(i);

  // for(i = 0; i < 10; ++i)
  //   njm::runif01();

  // std::vector<int> train,test;
  // std::priority_queue<std::pair<double,int> > queue;
  // for(i = 0; i < s.fD.numNodes; ++i)
  //   queue.push(std::pair<double,int>(njm::runif01(),i));

  // std::pair<double,int> top;
  // for(i = 0; i < s.fD.numNodes; ++i){
  //   top = queue.top();
  //   queue.pop();
    
  //   if(i < int(double(s.fD.numNodes)*0.6 + 1))
  //     train.push_back(top.second);
  //   else
  //     test.push_back(top.second);
  // }

  // std::cout << njm::toString(train," ","\n")
  // 	    << njm::toString(test," ","\n");


  // oq.qEval.buildRD(train);
  // oq.qEval.solve();
  // std::cout << oq.qEval.beta.sum() << std::endl;
  // std::cout << oq.qEval.bellRes() << std::endl;
  // oq.qEval.buildRD(test);
  // std::cout << oq.qEval.bellRes() << std::endl;
  // oq.qEval.buildRD(train);
  // std::cout << oq.qEval.bellRes() << std::endl;
  // oq.qEval.buildRD(test);
  // std::cout << oq.qEval.bellRes() << std::endl;
  
  njm::sett.clean();
  return 0;
}
