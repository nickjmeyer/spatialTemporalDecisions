#include "timer.hpp"

Timer njm::timer;

using namespace std;
using namespace chrono;

Timer::Timer(){
#ifndef NJM_NO_TIMER
  // initialize containers
  running.resize(omp_get_max_threads());
  tick.resize(omp_get_max_threads());
#endif
}


Timer::~Timer(){
#ifndef NJM_NO_TIMER
  // print the times before destroying
  print();
#endif
}


void Timer::print(){
#ifndef NJM_NO_TIMER
  // print the times
  
  int i,I = omp_get_max_threads();
  map<string,high_resolution_clock::duration> total;
  map<string,int> counts;

  // count up the times across all threads
  map<string,high_resolution_clock::duration>::iterator it,beg,end;
  for(i = 0; i < I; ++i){
    // iterate over threads
    
    beg = running.at(i).begin();
    end = running.at(i).end();
    for(it = beg; it != end; ++it){
      // iterate over chunks that the thread timed
      
      if(counts.find(it->first) == counts.end()){
	counts.insert(pair<string,int>(it->first,0));
	total.insert(pair<string,high_resolution_clock::duration>
		     (it->first,high_resolution_clock::duration::zero()));
      }

      // count how many threads contributed
      ++counts.at(it->first);

      // update total time
      total.at(it->first) += it->second;
    }
  }


  // sort the times
  priority_queue<pair<high_resolution_clock::duration,string> > totalOrd;
  beg = total.begin();
  end = total.end();

  high_resolution_clock::duration totalTime;
  totalTime = high_resolution_clock::duration::zero();
  for(it = beg; it != end; ++it){
    totalTime += it->second;
    totalOrd.push(pair<high_resolution_clock::duration,string>
		  (it->second,it->first));
  }


  // if at least one section print the times
  if(totalOrd.size() > 0){
    // header
    printf("%s%s","\n#-----------------------------------TIMER",
	   "------------------------------------\n");
    printf("%18s  %12s %18s %18s\n","ID","Percent","Total","Total/Thread");

    // print out the [name, percentage, total, total per thread]
    I = totalOrd.size();
    pair<high_resolution_clock::duration,string> top;
    double percTotal,totalMins,perThreadMins;
    milliseconds totalTimeMs,curTimeMs;
    totalTimeMs = duration_cast<milliseconds>(totalTime);
    for(i = 0; i < I; ++i){
      top = totalOrd.top();
      totalOrd.pop();

      curTimeMs = duration_cast<milliseconds>(top.first);

      percTotal = 100.0*double(curTimeMs.count()) / double(totalTimeMs.count());
      totalMins = double(curTimeMs.count()) * double(1.66666666667e-5);
      perThreadMins = totalMins / double(counts.at(top.second));
    
      printf("%18.18s: % 12.8f % 18.8f % 18.8f   mins\n",
	     top.second.c_str(), percTotal, totalMins, perThreadMins);
    }
    // footer
    printf("\n%18.18s: % 12.8f   mins\n",
	   "Aggregated total",
	   double(totalTimeMs.count()) * double(1.66666666667e-5));
    printf("%s%s","------------------------------------TIMER",
	   "-----------------------------------#\n");
  }
  
#endif
}


void Timer::start(const string name){
#ifndef NJM_NO_TIMER
  // start the time (dont get time till end of function)
  
  int thread = omp_get_thread_num();

  map<string,time_point<high_resolution_clock> >::iterator it,end;

  end = tick.at(thread).end();
  it = tick.at(thread).find(name);

  if(end == it){ // if this chunk doesn't have a time yet
    running.at(thread).insert(pair<string,high_resolution_clock::duration>
			      (name,high_resolution_clock::duration::zero()));
    tick.at(thread).insert(pair<string,time_point<high_resolution_clock> >
			   (name,
			    high_resolution_clock::now()));
  }
  else // otherwise overwrite the current tick
    it->second = high_resolution_clock::now();
#endif
}


void Timer::stop(const string name){
#ifndef NJM_NO_TIMER
  // get the time first thing
  time_point<high_resolution_clock> tock = high_resolution_clock::now();

  int thread = omp_get_thread_num();

  // add the duration
  high_resolution_clock::duration diff;
  diff = tock.time_since_epoch() - tick.at(thread).at(name).time_since_epoch();
  running.at(thread).at(name) += diff;
#endif
}
