#include "modelMulti.hpp"

double
MultiModel::oneOnOne(const int notNode,
		     const int infNode,
		     const SimData & sD,
		     const TrtData & tD,
		     const FixedData & fD,
		     const DynamicData & dD) const {
  return m.at(ind).oneOnOne(notNode,infNode,sD,tD,fD,dD);
}



void
MultiModel::fit(const SimData & sD, const TrtData & tD,
		const FixedData & fD, const DynamicData & dD,
		MultiParam & mP){
  int i,s = std::tuple_size(decltype(mods));
  for(i = 0; i < s; ++i){
    m.at(i).fitType = fitType;
    m.at(i).fit(sD,tD,fD,dD,mP.pars.at(i));
  }
}



void MultiModel::fit(const SimData & sD, const TrtData & tD,
		     const FixedData & fD,
		     const DynamicData & dD,
		     MultiParam & mP,
		     const MultiParam mPInit){
  int i,s = std::tuple_size(decltype(mods));
  for(i = 0; i < s; ++i){
    mods.at(i).fitType = fitType;
    mods.at(i).fit(sD,tD,fD,dD,mP.pars.at(i),mPInit.pars.at(i));
  }
}
