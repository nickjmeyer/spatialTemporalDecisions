#ifndef RUN_STATS_HPP
#define RUN_STATS_HPP

#include <vector>
#include <cmath>
#include <algorithm>

class RunStats {
protected:
    std::vector<double> vals;
    unsigned int n;

    std::vector<double> reservoir;
    unsigned int reservoir_cnt;

    double sMean_; // sample mean
    double sSqMean_; // sample mean of squared values
    double sVar_; // sample variance
    double sSd_; // sample standard deviation
    double seMean_; // standard error of mean

public:
    RunStats();
    RunStats(const std::vector<double> & init);
    RunStats(const unsigned int & num_vals);
    RunStats(const unsigned int & num_vals,
            const std::vector<double> & init);

    void update(const double & add);

    void update(const unsigned int & position,
            const double & add);

    void update(const std::vector<double> & add);

    void update(const std::vector<unsigned int> & positions,
            const std::vector<double> & add);

    void squash();

    double sMean() const;

    double sVar() const;

    double sSd() const;

    double seMean() const;
};




#endif
