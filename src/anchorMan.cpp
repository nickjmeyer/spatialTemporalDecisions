#include "anchorMan.hpp"


std::vector<double> AnchorManTunePar::getPar() const{
  return std::vector<double>();
}

void AnchorManTunePar(const std::vector<double> & par){
  // do nothing;
}


AnchorManTunePar::AnchorManTunePar(){
  numSamples = 100;
  cutoff = .95;
}
