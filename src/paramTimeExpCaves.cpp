#include "paramTimeExpCaves.hpp"


unsigned int ParamTimeExpCaves::initParsSize(const FixedData & fD){
  return 1U;
}


void ParamTimeExpCaves::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  time = std::vector<int>(numNodes,0);
  iPropCaves = fD.propCaves;
  std::for_each(iPropCaves.begin(),iPropCaves.end(),
		[](double & x){ x = 1.0/x; });
}


void ParamTimeExpCaves::updateBefore(){
}


void ParamTimeExpCaves::updateAfter(){
}


void ParamTimeExpCaves::setFill(std::vector<double> & probs,
				const SimData & sD,
				const TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD){
  double xi = *beg;
  
  int i,j;
  std::vector<double>::iterator it0;
  std::vector<double>::const_iterator it1, beg;

  std::vector<double> xiTime(numNodes,0);
  
  time = sD.timeInf;
  for(i = 0; i < numNodes; ++i){
    xiTime.at(i) = xi * std::exp(iPropCaves.at(i)*double(time.at(i)-1));
  }
  
  beg = xiTime.begin();
  for(i = 0, it0 = probs.begin(); i < numNodes; ++i){
    for(j = 0, it1 = beg; j < numNodes; ++j,++it0,++it1){
      *it0 += *it1;
    }
  }
}

void ParamTimeExpCaves::modFill(std::vector<double> & probs,
				const SimData & sD,
				const TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD){
  double xi = *beg;
  
  int i,j;
  std::vector<double>::iterator it0;

  std::vector<int> change;
  std::vector<double> diff;

  change.reserve(numNodes);
  diff.reserve(numNodes);

  int lapse;
  double val0, val1;
  
  for(i = 0; i < numNodes; ++i){
    lapse = sD.timeInf.at(i) - time.at(i);
    if(lapse != 0){
      change.push_back(1);
      val0 = std::exp(iPropCaves.at(i)*double(time.at(i)-1));
      val1 = std::exp(iPropCaves.at(i)*double(sD.timeInf.at(i)-1));
      diff.push_back(xi * (val1 - val0));
    }
    else{
      change.push_back(0);
      diff.push_back(0);
    }
  }
 
  for(i = 0, it0 = probs.begin(); i < numNodes; ++i){
    for(j = 0; j < numNodes; ++j,++it0){
      if(change.at(j) == 1){
	*it0 += diff.at(j);
      }
    }
  }

  time = sD.timeInf;
}
