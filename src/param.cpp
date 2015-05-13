#include "param.hpp"

ParamBase::ParamBase(){
}

ParamBase::ParamBase(const ParamBase & p){
  pars = p.pars;
  beg = pars.begin();
  end = pars.end();
  parsSize = p.parsSize;
  names = p.names;
}


void ParamBase::init(const FixedData & fD){
  // initialize the size
  parsSize = initParsSize(fD);

  // initialize vector to zero
  pars = std::vector<double> (parsSize,0);
  beg = pars.begin();
  end = pars.end();

  // names
  names = initNames();

  // initialize the internal information
  initInternal(fD);
}


void ParamBase::save(const std::string & m) const {
  unsigned int i;
  for(i = 0; i < parsSize; ++i)
    njm::toFile(pars.at(i),
		njm::sett.srcExt("./Param"+m+"/"+names.at(i)+".txt"),
		std::ios_base::out);
}


void ParamBase::read(const std::string & m){
  unsigned int i;
  double val;
  for(i = 0; i < parsSize; ++i){
    njm::fromFile(val,
		  njm::sett.srcExt("./Param"+m+"/"+names.at(i)+".txt"));
    pars.at(i) = val;
  }
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
  
