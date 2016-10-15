#include <iostream>
#include <list>
#include <vector>
#include <queue>
#include <cmath>
#include <utility>
#include <set>
#include <Rcpp.h>


// [[Rcpp::export]]
std::vector<int> getClusters(const std::vector<int> & net,
        const int numNodes);

// [[Rcpp::export]]
std::vector<int> mergeOneCluster(std::vector<int> net,
        const std::vector<int> & clusters,
        const std::vector<double> & x,
        const std::vector<double> & y,
        const int numNodes);


std::vector<int> getClusters(const std::vector<int> & net,
        const int numNodes){
    std::vector<int> clusters(numNodes,-1);

    std::list<int> left;
    int i;
    for(i = 0; i < numNodes; ++i){
        left.push_back(i);
    }

    int j;
    int curGroup = 0;
    while(!left.empty()){
        int node = left.front();

        std::vector<int> conn0(numNodes,0),conn1(numNodes,0);
        int cnt0 = 0, cnt1 = 0;
        for(i = 0; i < numNodes; ++i){
            if(net.at(node*numNodes + i) == 1 || i == node){
                conn0.at(i) = 1;
                ++cnt0;
            }
        }

        do{
            cnt1 = cnt0;
            for(i = 0; i < numNodes; ++i){
                if(conn0.at(i) && !conn1.at(i)){
                    conn1.at(i) = 1;
                    for(j = 0; j < numNodes; ++j){
                        if(net.at(i*numNodes + j) && !conn0.at(j)){
                            conn0.at(j) = 1;
                            ++cnt0;
                        }
                    }
                }
            }
        }while(cnt1!=cnt0);

        for(i = 0; i < numNodes; ++i){
            if(conn0.at(i)){
                clusters.at(i) = curGroup;
                left.remove(i);
            }
        }

        ++curGroup;
    }

    return clusters;
}



std::vector<int> mergeOneCluster(std::vector<int> net,
        const std::vector<int> & clusters,
        const std::vector<double> & x,
        const std::vector<double> & y,
        const int numNodes){
    std::set<int> names;
    int i;
    for(i = 0; i < numNodes; ++i)
        names.insert(clusters.at(i));

    std::vector<std::pair<int,int> > combos;

    std::set<int>::const_iterator it0,it1,beg,end;
    beg = names.begin();
    end = names.end();
    for(it0 = beg; it0 != end; ++it0){
        for(it1 = it0; it1 != end; ++it1){
            if(it1 != it0)
                combos.push_back(std::pair<int,int>(*it0,*it1));
        }
    }



    std::priority_queue<std::pair<double,std::pair<int,int> > > pairs;

    int numCombos = combos.size();
    int j,k;
    for(i = 0; i < numCombos; ++i){
        int c0,c1;
        c0 = combos.at(i).first;
        c1 = combos.at(i).second;

        for(j = 0; j < numNodes; ++j){
            if(clusters.at(j) == c0){
                for(k = 0; k < numNodes; ++k){
                    if(clusters.at(k) == c1){
                        double dist = std::sqrt(std::pow(x.at(j) - x.at(k),2.0) +
                                std::pow(y.at(j) - y.at(k),2.0));
                        std::pair<int,int> nodes(j,k);
                        pairs.push(std::pair<double,std::pair<int,int> >(-dist,nodes));
                    }
                }
            }
        }
    }

    std::pair<int,int> nodes = pairs.top().second;
    net.at(nodes.first*numNodes + nodes.second) = 1;
    net.at(nodes.second*numNodes + nodes.first) = 1;

    return(net);
}
