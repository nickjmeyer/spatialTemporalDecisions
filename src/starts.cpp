#include "starts.hpp"

Starts::Starts(const int numReps, const int numNodes, const int dynamic){
  int num = std::max(numNodes/100,1);
  this->dynamic = dynamic;

  int i,j;
  ind.resize(numReps);
  std::pair<double,int> top;
  for(i = 0; i < numReps; ++i){
    std::priority_queue<std::pair<double,int> > q;
    for(j = 0; j < numNodes; ++j){
      q.push(std::pair<double,int>(njm::runif01(),j));
    }

    ind.at(i).clear();
    for(j = 0; j < num; ++j){
      top = q.top();
      q.pop();

      ind.at(i).push_back(top.second);
    }
  }
}


std::vector<int> Starts::operator[](const int i){
  if(dynamic)
    return ind.at(i);
  else
    return ind.at(0);
}
