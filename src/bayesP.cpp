#include "bayesP.hpp"

std::vector<double> getStats(const std::vector<std::vector<int> > & h,
			     const SimData & sD,
			     const FixedData & fD){
  std::vector<double> stats;
  // total infected
  stats.push_back(sD.numInfected);

  // number of newly infected by year
  // int i,j,sum = 0,pastSum = 0;
  // for(i = 0; i < (int)h.size(); ++i){
  //   pastSum = sum;
  //   sum = 0;
  //   for(j = 0; j < fD.numNodes; ++j)
  //     if(h.at(i).at(j) >= 2)
  // 	++sum;

  //   if(i > 0)
  //     stats.push_back(sum - pastSum);
  // }
  int i,j,sum = 0;
  for(i = 1; i < (int)h.size(); ++i){
    sum = 0;
    for(j = 0; j < fD.numNodes; ++j)
      if(h.at(i).at(j) >= 2 && h.at(i-1).at(j) < 2)
  	++sum;

    stats.push_back(sum);
  }


  // average year of infection
  sum = 0;
  for(i = 0; i < (int)h.size(); ++i){
    for(j = 0; j < fD.numNodes; ++j){
      if((i == 0 && h.at(i).at(j) >= 2) ||
	 (i > 0 && h.at(i).at(j) >= 2 && h.at(i-1).at(j) < 2))
	sum += i;
    }
  }
  stats.push_back(((double)sum)/((double)sD.numInfected));


  // average lat long for infected
  double cLong,cLat;
  cLong = cLat = 0;
  for(i = 0 ; i < sD.numInfected; ++i){
    cLong += fD.centroidsLong.at(sD.infected.at(i));
    cLat += fD.centroidsLat.at(sD.infected.at(i));
  }
  stats.push_back(cLong/((double)sD.numInfected));
  stats.push_back(cLat/((double)sD.numInfected));

  // average spread distance from start for infected
  double sDist;
  int start = -1,numStart = 0;
  for(i = 0; i < fD.numNodes; ++i){
    if(h.at(0).at(i) >= 2){
      if(numStart == 0)
	start = i;
      ++numStart;
    }
  }
  if(numStart > 1)
    std::cout << "warning...multiple starting locations are listed...."
	      << std::endl;
  sDist = 0;
  for(i = 0; i < sD.numInfected; ++i){
    sDist += fD.dist.at(sD.infected.at(i)*fD.numNodes + start);
  }
  sDist/=(double)sD.numInfected;
  stats.push_back(sDist);

  // min lat long for infected
  cLong = std::numeric_limits<double>::max();
  cLat = std::numeric_limits<double>::max();
  for(i = 0; i < sD.numInfected; ++i){
    if(fD.centroidsLong.at(sD.infected.at(i)) < cLong)
      cLong = fD.centroidsLong.at(sD.infected.at(i));
    if(fD.centroidsLat.at(sD.infected.at(i)) < cLat)
      cLat = fD.centroidsLat.at(sD.infected.at(i));
  }
  stats.push_back(cLong);
  stats.push_back(cLat);


  // max lat long for infected
  cLong = std::numeric_limits<double>::lowest();
  cLat = std::numeric_limits<double>::lowest();
  for(i = 0; i < sD.numInfected; ++i){
    if(fD.centroidsLong.at(sD.infected.at(i)) > cLong)
      cLong = fD.centroidsLong.at(sD.infected.at(i));
    if(fD.centroidsLat.at(sD.infected.at(i)) > cLat)
      cLat = fD.centroidsLat.at(sD.infected.at(i));
  }
  stats.push_back(cLong);
  stats.push_back(cLat);

  
  // max dist from starting location
  sDist = std::numeric_limits<double>::lowest();
  for(i = 0; i < sD.numInfected; ++i)
    if(fD.dist.at(sD.infected.at(i)*fD.numNodes + start) > sDist)
      sDist = fD.dist.at(sD.infected.at(i)*fD.numNodes + start);
  stats.push_back(sDist);


  return stats;
}



int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityTimeInfModel GM;
  typedef GravityTimeInfParam GP;
  typedef GM EM;
  typedef GP EP;

  typedef System<GM,GP,EM,EP> S;

  S s("obsData.txt");

  std::vector<std::vector<int> > h;
  h = s.sD.history;
  h.push_back(s.sD.status);


  njm::message(getStats(h,s.sD,s.fD));

  return 0;
}
