#ifndef TIMER_HPP__
#define TIMER_HPP__


#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <string>
#include <chrono>
#include <omp.h>


class Timer {
 public:
  Timer();
  ~Timer();

  void start(const std::string name);
  void stop(const std::string name);

 private:
  std::vector<std::map<std::string,std::chrono::milliseconds> > running;
  std::vector<std::map<std::string,
		       std::chrono::time_point<
			 std::chrono::high_resolution_clock> > > tick;
};

namespace njm {
	extern Timer timer;
}
  


#endif
