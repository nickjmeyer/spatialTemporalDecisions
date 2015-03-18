#include "timer.hpp"

Timer njm::timer;

using namespace std;
using namespace chrono;

Timer::Timer(){
#ifndef NJM_NO_TIMER
  running.resize(omp_get_max_threads());
  tick.resize(omp_get_max_threads());
#endif
}


Timer::~Timer(){
#ifndef NJM_NO_TIMER

  int i,I = omp_get_max_threads();
  map<string,milliseconds> total;
  map<string,int> counts;

  map<string,milliseconds>::iterator it,beg,end;
  for(i = 0; i < I; ++i){
    beg = running.at(i).begin();
    end = running.at(i).end();
    for(it = beg; it != end; ++it){
      if(counts.find(it->first) == counts.end()){
	counts.insert(pair<string,int>(it->first,0));
	total.insert(pair<string,milliseconds>(it->first,milliseconds::zero()));
      }

      ++counts.at(it->first);

      total.at(it->first) += it->second;
    }
  }


  priority_queue<pair<milliseconds,string> > totalOrd;
  beg = total.begin();
  end = total.end();

  milliseconds totalTime = milliseconds::zero();
  for(it = beg; it != end; ++it){
    totalTime += it->second;
    totalOrd.push(pair<milliseconds,string>(it->second,it->first));
  }


  printf("%s%s","\n#-----------------------------------TIMER",
	 "------------------------------------\n");
  printf("%18s  %12s %18s %18s\n","ID","Percent","Total","Total/Thread");

  I = totalOrd.size();
  pair<milliseconds,string> top;
  double percTotal,totalMins,perThreadMins;
  for(i = 0; i < I; ++i){
    top = totalOrd.top();
    totalOrd.pop();

    percTotal = (double(top.first.count()) / double(totalTime.count()))*100.0;
    totalMins = double(top.first.count()) * double(1.66666666667e-5);
    perThreadMins = totalMins / double(counts.at(top.second));
    
    printf("%18.18s: % 12.8f % 18.8f % 18.8f   mins\n",
	   top.second.c_str(), percTotal, totalMins, perThreadMins);
  }
  printf("%s%s","------------------------------------TIMER",
	 "-----------------------------------#\n");
  
#endif
}


void Timer::start(const string name){
#ifndef NJM_NO_TIMER
  int thread = omp_get_thread_num();

  map<string,time_point<high_resolution_clock> >::iterator it,end;

  end = tick.at(thread).end();
  it = tick.at(thread).find(name);

  if(end == it){
    running.at(thread).insert(pair<string,milliseconds>
			      (name,milliseconds::zero()));
    tick.at(thread).insert(pair<string,time_point<high_resolution_clock> >
			   (name, high_resolution_clock::now()));
  }
  else
    it->second = high_resolution_clock::now();
#endif
}


void Timer::stop(const string name){
#ifndef NJM_NO_TIMER
  time_point<high_resolution_clock> tock = high_resolution_clock::now();

  int thread = omp_get_thread_num();
  
  milliseconds diff;
  diff = duration_cast<milliseconds>
    (tock.time_since_epoch() - tick.at(thread).at(name).time_since_epoch());

  running.at(thread).at(name) += diff;
#endif
}
