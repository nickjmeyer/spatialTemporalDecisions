#include "modelMulti.hpp"


MultiModel::MultiModel(){
  fill();

  if(int(m.size()) != numModels){
    std::cout << "MultiModel: number of models is not correct."
	      << std::endl;
    throw(1);
  }
}


MultiModel::MultiModel(const MultiModel & mm){
  fill();
  ind = mm.ind;
  setType(mm.getType());
  int s, S = size();
  for(s = 0; s < S; ++s)
    (*m.at(s)) = (*mm.m.at(s));
  
  if(int(m.size()) != numModels){
    std::cout << "MultiModel: number of models is not correct."
	      << std::endl;
    throw(1);
  }
}


MultiModel::~MultiModel(){
  int i,s = m.size();
  for(i = 0; i < s; ++i){
    delete m.at(i);
  }
}


MultiModel & MultiModel::operator=(const MultiModel & mm){
  if(this != & mm){
    this->MultiModel::~MultiModel();
    new (this) MultiModel(mm);
  }
  return *this;
}


void MultiModel::fill(){
  // ***************************************** //
  // make sure to change numModels in hpp file
  // ***************************************** //
  m.clear();
  m.push_back(new RadiusModel);
  // m.push_back(new CaveModel);
  // m.push_back(new RangeModel);
  // m.push_back(new GravityTimeInfExpCavesModel);
  m.push_back(new GravityModel);
}


void MultiModel::setType(const Estimation & est){
  getType() = est;
  int s,S = size();
  for(s = 0; s < S; ++s)
    m.at(s)->setType(est);
}


void MultiModel::load(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD){
  int i,s = size();
  for(i = 0; i < s; ++i)
    m.at(i)->load(sD,tD,fD,dD);
}

void MultiModel::infProbs(const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD){
  int i,s = size();
  for(i = 0; i < s; ++i)
    m.at(i)->infProbs(sD,tD,fD,dD);
}
  
void MultiModel::update(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD){
  int i,s = size();
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
  int i,s = size();
  for(i = 0; i < s; ++i){
    m.at(i)->fit(sD,tD,fD,dD,useInit);
  }
}


int MultiModel::size(){
  return numModels;
}


void MultiModel::modSel(const int & ind){
  this->ind = ind;
}


BaseModel * MultiModel::operator[](const int i){
  return m.at(i);
}
