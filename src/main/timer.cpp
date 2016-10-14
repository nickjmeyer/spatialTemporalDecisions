#include "timer.hpp"

namespace njm{

Timer timer;

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
    map<string,DurAndCalls> total;

    // count up the times across all threads
    map<string,DurAndCalls>::iterator it,beg,end;
    for(i = 0; i < I; ++i){
        // iterate over threads

        beg = running.at(i).begin();
        end = running.at(i).end();
        for(it = beg; it != end; ++it){
            // iterate over chunks that the thread timed

            if(total.find(it->first) == total.end()){
                total.insert(pair<string,DurAndCalls>(
                                it->first,DurAndCalls(
                                        high_resolution_clock::duration::zero(),
                                        0)));
            }

            // update total time
            total.at(it->first).first += it->second.first;
            // update total calls
            total.at(it->first).second += it->second.second;
        }
    }


    // sort the times
    priority_queue<pair<DurAndCalls, string> > totalOrd;
    beg = total.begin();
    end = total.end();

    high_resolution_clock::duration totalTime;
    totalTime = high_resolution_clock::duration::zero();
    for(it = beg; it != end; ++it){
        totalTime += it->second.first;
        totalOrd.push(pair<DurAndCalls, string>(
                        DurAndCalls(it->second.first,it->second.second),
                        it->first));
    }


    // if at least one section print the times
    if(totalOrd.size() > 0){
        // header
        printf("%s%s","\n#-----------------------------------TIMER",
                "------------------------------------\n");
        printf("%18s  %12s %18s %18s\n","ID","Calls","Total","Total/Call");

        // print out the [name, calls, total, total per call]
        I = totalOrd.size();
        pair<DurAndCalls,string> top;
        double totalMins,perCallMins;
        milliseconds totalTimeMs,curTimeMs;
        unsigned long long calls;
        totalTimeMs = duration_cast<milliseconds>(totalTime);
        for(i = 0; i < I; ++i){
            top = totalOrd.top();
            totalOrd.pop();

            curTimeMs = duration_cast<milliseconds>(top.first.first);
            totalMins = double(curTimeMs.count()) * double(1.66666666667e-5);
            calls = top.first.second;
            perCallMins = totalMins / double(calls);

            printf("%18.18s: %012llu % 18.8f % 18.8f   mins\n",
                    top.second.c_str(), calls, totalMins, perCallMins);
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
        running.at(thread).insert(pair<string,DurAndCalls>
                (name,DurAndCalls(high_resolution_clock::duration::zero(),0)));
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
    running.at(thread).at(name).first += diff;
    ++running.at(thread).at(name).second;
#endif
}

}
