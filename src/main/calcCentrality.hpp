#ifndef CALC_CENTRALITY_HPP
#define CALC_CENTRALITY_HPP

#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <armadillo>
// #include <eigen3/Eigen/Sparse>

void getSubGraph(int nodes, std::vector<int> const * const network,
        std::vector<double> * const subGraph);
void getSubGraph(int nodes, std::vector<int> const * const network,
        std::vector<double> * const subGraph, int const deg);
// void getSubGraphEig(int nodes, std::vector<int> const * const network,
// 		    std::vector<double> * const subGraph, int const deg);
void getBetweenness(int nodes, std::vector<int> const * const network,
        std::vector<double> * const btwn);

class MyQueue{
public:
    int pop();
    void push(int v);
    void clear();
    int empty();
private:
    std::list<int> vals;
};

class MyStack{
public:
    int pop();
    void push(int v);
    void clear();
    int empty();
private:
    std::list<int> vals;
};


#endif
