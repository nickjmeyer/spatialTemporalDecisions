#ifndef TIMER_HPP
#define TIMER_HPP


#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <string>
#include <chrono>
#include <omp.h>

using namespace std::chrono;

typedef std::pair<high_resolution_clock::duration,
                  unsigned long long>
DurAndCalls;

namespace njm {

class Timer {
public:
    Timer();
    ~Timer();

    void print();

    void start(const std::string name);
    void stop(const std::string name);

private:
    std::vector<std::map<std::string, DurAndCalls> > running;

    std::vector<std::map<std::string,
                         time_point<high_resolution_clock,
                                    high_resolution_clock::duration> > > tick;
};

extern Timer timer;
}



#endif
