#include "starts.hpp"

Starts::Starts(const std::string & file){
    dynamic = 0;

    std::vector<int> start;
    njm::fromFile(start,njm::sett.srcExt(file));
    ind.clear();
    ind.push_back(start);
}

Starts::Starts(const std::vector<int> & start) {
    dynamic = 0;

    ind.clear();
    ind.push_back(start);
}

Starts::Starts(const int numReps, const int numNodes){
    njm::resetSeed(0);

    dynamic = 1;

    int num = std::max(numNodes/100,1);

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

        // increasing order
        std::sort(ind.at(i).begin(),ind.at(i).end());
    }
}


std::vector<int> Starts::operator[](const int i) const{
    if(dynamic)
        return ind.at(i);
    else
        return ind.at(0);
}
