#include "test.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef GravityTimeInfExpCavesModel GM;
    typedef GravityTimeInfExpCavesParam GP;
    typedef GM EM;
    typedef GP EP;

    typedef System<GM,GP,EM,EP> S;

    typedef ToyFeatures2<EM,EP> F;
    typedef RankAgent<F,EM,EP> RA;

    S s;
    RA ra;

    ra.reset();

    s.modelGen.fitType = MLE;
    s.modelEst.fitType = MLE;

    int i;

    for(i = 0 ; i < 10; ++i)
      njm::runif01();
    
    for(i = 0 ; i < (s.fD.finalT+10); ++i){
      if(i >= s.fD.trtStart)
      	ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelGen,s.paramGen);
      s.nextPoint();
    }

    std::cout << s.value() << std::endl;
    
  }


  njm::sett.clean();
  return 0;
}
