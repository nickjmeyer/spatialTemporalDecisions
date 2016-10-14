#ifndef SORT_MERGE
#define SORT_MERGE

#include <vector>
#include <map>
#include "rand.hpp"

template<typename T>
void sortByValue(std::vector<T> * const A);
// REQUIRES: nothing
// MODIFIES: *A
// EFFECTS: sorts A in ascending order


template<typename T>
void mergeByValue(std::vector<T> const * const A,
		  std::vector<T> * const B,
		  int const iLeft, int const iRight, int const iEnd);
// REQUIRES: nothing
// MODIFIES: *B
// EFFECTS: merges sections of *A into *B by comparing values of *A


void sortByVec(std::vector<int> * const A,
	       std::vector<double> * const int2double);
// REQUIRES: all elements of *A are keys in int2double
// MODIFIES: *A
// EFFECTS: sorts A in ascending order by comparing elements in int2double
//          using sort-merge technique.  This calls sortByMap




void sortByMap(std::vector<int> * const A,
	       std::map<int, double> const * const int2double);
// REQUIRES: all elements of *A are keys in int2double
// MODIFIES: *A
// EFFECTS: sorts A in ascending order by comparing keys in int2double using
//          the sort-merge technique

void mergeByMap(std::vector<int> const * const A,
		std::map<int,double> const * const  int2double,
		std::vector<int> * const B,
		int const iLeft, int const iRight, int const iEnd);
// REQUIRES: all elements of *A are in int2double
// MODIFIES: *B
// EFFECTS: merges sections of *A into *B by comparing keys in int2double


#endif
