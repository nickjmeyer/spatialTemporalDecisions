#include "tuneSp.hpp"


FFX::FFX(){
  numFactor = 0;
  numCombo = 1;
  numReps = 10;
}


void FFX::addFactor(const std::string & f,
		    const std::vector<double> & fVals){
  factors.push_back(f);
  values.push_back(fVals);
  maxSett.push_back((int)fVals.size());

  numCombo *= (int)fVals.size();
  ++numFactor;
}

void FFX::addStat(const std::string & s){
  stats.push_back(s);
}


std::vector<double> FFX::getSett(const int i) const {
  if(i > (numCombo * numReps)){
    std::cout << "invalid setting index" << std::endl;
    throw(1);
  }

  int i,maxInd=numCombo,ind = i % numCombo;
  std::vector<double> sett;
  for(i = 0 ; i < numFactor; i++){
    maxInd /= (int)maxSett.at(i);
    sett.push_back(values.at(i).at(ind/maxInd));
    ind = ind % maxInd;
    ind /= maxInd;

    
  }
}


int main(){


  std::vector<double> Avals = {10,30,50};
  std::vector<double> Bvals = {0.1,0.3,0.5};
  std::vector<double> Cvals = {2.0,5.0,10.0};
  std::vector<double> Tvals = {0.2,0.5,1.0};
  std::vector<double> Lvals = {0.5,1.0,2.0};


  
  return 0;
}

