#include "anchorMan.hpp"


template <class System, class Agent>
AnchorMan::AnchorMan(){
}



template <class System, class Agent>
void AnchorMan::addPar(const std::vector<double> & par){
  parHist.push_back(par);
}



