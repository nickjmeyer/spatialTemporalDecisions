#ifndef STARTS_HPP
#define STARTS_HPP


#include <vector>
#include <queue>
#include "utilities.hpp"
#include "settings.hpp"
#include "rand.hpp"


class Starts {
public:
    Starts(const std::string & file);
    Starts(const int numReps, const int numNodes);

    std::vector<int> operator[](const int i) const;

private:
    int dynamic;
    std::vector<std::vector<int> > ind;
};


#endif
