#include "param.hpp"

ParamBase::ParamBase(const ParamBase & p){
  pars = p.pars;
  beg = pars.begin();
  end = pars.end();
  parsSize = p.parsSize;
}


void ParamBase::init(const FixedData & fD){
  // scalar parameter

  // initialize the size
  parsSize = initParsSize(fD);

  // initialize vector to zero
  pars = std::vector<double> (parsSize,0);
  beg = pars.begin();
  end = pars.end();

  // initialize the internal information
  initInternal(fD);
}


std::vector<double> ParamBase::getPar() const {
  return pars;
}


std::vector<double>::const_iterator
ParamBase::putPar(std::vector<double>::const_iterator newParIt){
  updateBefore();
  std::vector<double>::iterator it;
  for(it = beg; it != end; ++it, ++newParIt)
    *it = *newParIt;
  updateAfter();
  return newParIt;
}


unsigned int ParamBase::size() const {
  return parsSize;
}


std::vector<double> ParamBase::partial2(const int notNode,
					const int infNode,
					const SimData & sD,
					const TrtData & tD,
					const FixedData & fD,
					const DynamicData & dD){
  return std::vector<double>(parsSize*parsSize,0);
}
  
