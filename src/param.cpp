#include "param.hpp"


BaseParam::BaseParam(const FixedData & fD){
  // scalar parameter

  // initialize the size
  parsSize = initParsSize(fD);

  // initialize vector to zero
  std::vector<double> pars(parsSize,0);
  beg = pars.begin();
  end = pars.end();

  // initialize the internal information
  initInternal(fD);
}


std::vector<double> BaseParam::getPar() const {
  return pars;
}


void BaseParam::putPar(std::vector<double>::iterator & newParIt){
  std::vector<double>:::iterator it;
  for(it = beg; it != end; ++it, ++newParIt)
    *it = *newParIt;
  updateIternal();
}



