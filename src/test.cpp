#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCaves MG;
  typedef MG ME;
  typedef System<MG,ME> S;
  typedef ProximalAgent<ME> PA;

  S s;
  PA pa;

  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  s.modelEst_r = s.modelGen_r;
  s.revert();

  Starts starts("startingLocations.txt");
  s.reset(starts[0]);
  s.revert();

  std::vector<DataBundle> info;

  int i;
  for(i = 0; i < s.fD.finalT; ++i){
    if(i >= s.fD.trtStart)
      pa.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
    s.updateStatus();

    info.push_back(DataBundle(s.sD,s.tD,s.dD));

    s.nextPoint();
  }

  info.push_back(DataBundle(s.sD,s.tD,s.dD));

  std::vector<DataBundle> recov;

  std::vector<std::vector<int> > hist;
  hist = s.sD.history;
  hist.push_back(s.sD.status);
  
  recov = historyToData(hist);

  int equal = 1;
  for(i = 0; i < int(recov.size()) && equal == 1; ++i){
    std::cout << i << std::endl;
    if(equal == 1){
      equal = (std::get<0>(recov[i]) == std::get<0>(info[i]));
    }
    if(equal == 1){
      equal = (std::get<1>(recov[i]) == std::get<1>(info[i]));
    }
    if(equal == 1){
      equal = (std::get<2>(recov[i]) == std::get<2>(info[i]));
    }
  }
  std::cout << "equal: " << equal << std::endl;

  njm::sett.clean();
  return 0;
}
