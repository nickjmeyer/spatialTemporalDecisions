#include "paramTime.hpp"


unsigned int ParamTime::initParsSize(const FixedData & fD){
  return 1U;
}


void ParamTime::initInternal(const FixedData & fD){
  numNodes = fD.numNodes;
  time = std::vector<int>(numNodes,0);
}


void ParamTime::updateBefore(){
}


void ParamTime::updateAfter(){
}


void ParamTime::setFill(std::vector<double> & probs,
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
    xiTime.at(i) = xi * double(time.at(i) - 1);
  }
  
  beg = xiTime.begin();
  for(i = 0, it0 = probs.begin(); i < numNodes; ++i){
    for(j = 0, it1 = beg; j < numNodes; ++j,++it0,++it1){
      *it0 += *it1;
    }
  }
}

void ParamTime::modFill(std::vector<double> & probs,
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
  
  for(i = 0; i < numNodes; ++i){
    lapse = sD.timeInf.at(i) - time.at(i);
    if(lapse != 0){
      change.push_back(1);
      diff.push_back(xi * double(lapse));
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


std::vector<double> ParamTime::partial(const int notNode,
				       const int infNode,
				       const SimData & sD,
				       const TrtData & tD,
				       const FixedData & fD,
				       const DynamicData & dD){
  // use the information in the arguments since this will be called
  // for all past history information
  return std::vector<double>(1,double(sD.timeInf[infNode] - 1));
}
