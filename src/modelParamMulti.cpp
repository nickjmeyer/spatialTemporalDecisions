#include "modelParamMulti.hpp"


void MultiParam::load(){
  std::cout << "MultiParam: load() not defined" << std::endl;
  throw(1);
}

void MultiParam::save(){
  std::cout << "MultiParam: save() not defined" << std::endl;
  throw(1);
}


std::vector<double> MultiParam::getPar() const {
  std::cout << "MultiParam: getPar() not defined" << std::endl;
  throw(1);
  return std::vector<double>();
}

void MultiParam::putPar(const std::vector<double> & param){
  std::cout << "MultiParam: putPar() not defined" << std::endl;
  throw(1);
}




inline void MultiParam::setAll(){
  int i,s = std::tuple_size(decltype(pars));
  for(i = 0; i < s; ++i)
    std::get<i>(pars).setAll();
}


inline void MultiParam::setRow(const int r){
  int i,s = std::tuple_size(decltype(pars));
  for(i = 0; i < s; ++i)
    std::get<i>(pars).setRow(r);
}
inline void MultiParam::setCol(const int c){
  int i,s = std::tuple_size(decltype(pars));
  for(i = 0; i < s; ++i)
    std::get<i>(pars).setCol(c);
}
inline void MultiParam::setInd(const int r, const int c){
  int i,s = std::tuple_size(decltype(pars));
  for(i = 0; i < s; ++i)
    std::get<i>(pars).setInd(r,c);
}
