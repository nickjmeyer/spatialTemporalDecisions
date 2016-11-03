#include "runStats.hpp"
#include <glog/logging.h>

RunStats::RunStats()
    : n(0), sMean_(0.), sSqMean_(0.),
      sVar_(0.), sSd_(0.), seMean_(0.),
      reservoir_cnt(0) {
}

RunStats::RunStats(const std::vector<double> & init)
    : RunStats() {
    update(init);
}

void RunStats::update(const double & add){
    if(n == 0){
        ++n;
        sMean_ = add;
        sSqMean_ = add*add;
        sVar_ = 0.0;
        sSd_ = 0.0;
        seMean_ = 0.0;
    }
    else{
        ++n;
        sMean_ = sMean_ + (add - sMean_)/double(n);
        sSqMean_ = sSqMean_ + (add*add - sSqMean_)/double(n);
        sVar_ = (sSqMean_ - sMean_*sMean_)*double(n)/double(n-1);
        sSd_ = std::sqrt(sVar_);
        seMean_ = sSd_/std::sqrt(double(n));
    }
}


void RunStats::update(const unsigned int & position,
        const double & add){
    // make sure reservoir can fit a value at index position
    if (position >= reservoir.size()) {
        reservoir.resize(position+1);
    }

    reservoir.at(position) = add;
    ++reservoir_cnt;
}

void RunStats::update(const std::vector<double> & add){
    std::for_each(add.begin(),add.end(),
            [this](const double & x){
                update(x);
            });
}


void RunStats::update(const std::vector<unsigned int> & positions,
        const std::vector<double> & add){
    CHECK_EQ(positions.size(),add.size());
    const unsigned int size = positions.size();
    for (unsigned int i = 0; i < size; ++i) {
        update(positions.at(i),add.at(i));
    }
}


void RunStats::squash() {
    CHECK_EQ(reservoir_cnt,reservoir.size());

    std::vector<double>::const_iterator it,end;
    end = reservoir.end();
    for (it = reservoir.begin(); it != end; ++it) {
        update(*it);
    }
    reservoir.clear();
    reservoir_cnt = 0;
}


double RunStats::sMean() const {
    CHECK_EQ(reservoir_cnt,0)
        << "There are values in reservoir.  "
        << "Squash before calculating statistics.";
    return sMean_;
}


double RunStats::sVar() const {
    CHECK_EQ(reservoir_cnt,0)
        << "There are values in reservoir.  "
        << "Squash before calculating statistics.";
    return sVar_;
}


double RunStats::sSd() const {
    CHECK_EQ(reservoir_cnt,0)
        << "There are values in reservoir.  "
        << "Squash before calculating statistics.";
    return sSd_;
}


double RunStats::seMean() const {
    CHECK_EQ(reservoir_cnt,0)
        << "There are values in reservoir.  "
        << "Squash before calculating statistics.";
    return seMean_;
}
