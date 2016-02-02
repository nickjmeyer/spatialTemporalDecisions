#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef Model2GPowGDist M;
  typedef System<M,M> S;


  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  std::vector<std::string> names = {"power"};
  double power = s.modelGen_r.getPar(names)[0];
  s.modelGen_r.setPar("power",std::log(power));

  std::vector<double> par = s.modelGen_r.getPar();
  s.modelEst_r.putPar(par.begin());
  s.revert();

  Starts starts(15,s.fD.numNodes);
  s.reset(starts[3]);


  int i;
  for(i = 0; i < s.fD.trtStart; ++i){
    s.nextPoint();
  }

  std::cout << "numInfected: " << s.sD.numInfected << std::endl;

  s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,false);

  std::cout << "gen: " << njm::toString(s.modelGen.getPar()) << std::endl;
  std::cout << "est: " << njm::toString(s.modelEst.getPar()) << std::endl;

  for(i = 0; i < 10; i++){
    s.modelEst.sample();
    std::cout << i << ": " << njm::toString(s.modelEst.getPar()) << std::endl;
    std::cout << "==========================================" << std::endl;
  }

  njm::sett.clean();
  return 0;
}
