#include "modelMulti.hpp"


void MultiModel::load(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    m.at(i).load(sD,tD,fD,dD);
}

void MultiModel::infProbs(const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    m.at(i).infProbs(sD,tD,fD,dD);
}
  
void MultiModel::update(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    m.at(i).update(sD,tD,fD,dD);
}
  


double MultiModel::oneOnOne(const int notNode,
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
		const int & useInit){
  int i,s = m.size();
  for(i = 0; i < s; ++i){
    m.at(i).fitType = fitType;
    m.at(i).fit(sD,tD,fD,dD,useInit);
  }
}



void MultiModel::modSel(const int & ind){
  this->ind = ind;
}
