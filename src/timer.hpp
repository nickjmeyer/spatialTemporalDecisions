#ifndef TIMER_HPP__
#define TIMER_HPP__


#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <string>
#include <chrono>
#include <omp.h>

using namespace std::chrono;

class Timer {
 public:
  Timer();
  ~Timer();

  void print();
  
  void start(const std::string name);
  void stop(const std::string name);

 private:
  std::vector<std::map<std::string,high_resolution_clock::duration> > running;
  std::vector<std::map<std::string,
		       time_point<high_resolution_clock,
				  high_resolution_clock::duration> > > tick;
};

namespace njm {
	extern Timer timer;
}
  


#endif
