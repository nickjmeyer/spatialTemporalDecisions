#include "modelMulti.hpp"


MultiModel::MultiModel(){
  m.push_back(new GravityTimeInfExpModel);
  m.push_back(new GravityTimeInfExpCavesModel);
}

MultiModel::~MultiModel(){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    delete m.at(i);
}


void MultiModel::load(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    m.at(i)->load(sD,tD,fD,dD);
}

void MultiModel::infProbs(const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    m.at(i)->infProbs(sD,tD,fD,dD);
}
  
void MultiModel::update(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  int i,s = m.size();
  for(i = 0; i < s; ++i)
    m.at(i)->update(sD,tD,fD,dD);
}
  


double MultiModel::oneOnOne(const int notNode,
			    const int infNode,
			    const SimData & sD,
			    const TrtData & tD,
			    const FixedData & fD,
			    const DynamicData & dD) const {
  return m.at(ind)->oneOnOne(notNode,infNode,sD,tD,fD,dD);
}



void
MultiModel::fit(const SimData & sD, const TrtData & tD,
		const FixedData & fD, const DynamicData & dD,
		const int & useInit){
  int i,s = m.size();
  for(i = 0; i < s; ++i){
    m.at(i)->fitType = fitType;
    m.at(i)->fit(sD,tD,fD,dD,useInit);
  }
}


int MultiModel::size(){
  return int(m.size());
}


void MultiModel::modSel(const int & ind){
  this->ind = ind;
}


BaseModel * MultiModel::operator[](const int i){
  return m.at(i);
}
