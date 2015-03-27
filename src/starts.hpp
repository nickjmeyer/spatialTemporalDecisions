#ifndef STARTS_HPP__
#define STARTS_HPP__


#include <vector>
#include <queue>
#include "rand.hpp"


class Starts {
 public:
  Starts(const int numReps, const int numNodes, const int dynamic);

  std::vector<int> operator[](const int i);

 private:
  int dynamic;
  std::vector<std::vector<int> > ind;
};


#endif
