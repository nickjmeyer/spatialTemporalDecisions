#ifndef PARAM_HPP
#define PARAM_HPP

#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Eigen>
#include <string>
#include <boost/filesystem.hpp>
#include "data.hpp"


class ParamBase {
protected:

    std::vector<double> pars;
    std::vector<double>::iterator beg,end;
    unsigned int offset;
    unsigned int parsSize;

    std::vector<std::string> names;

    std::vector<bool> toScale;


    virtual unsigned int initParsSize(const FixedData & fD) = 0;

    virtual std::vector<std::string> initNames() = 0;

    virtual std::vector<bool> initToScale();

    // initialize the internal book-keeping data
    virtual void initInternal(const FixedData & fD) = 0;

    // update the internal book-keeping data
    // called inside of putPar() to guarantee up-to-date book-keeping data
    virtual void updateBefore() = 0;
    virtual void updateAfter() = 0;

public:
    ParamBase();
    ParamBase(const ParamBase & p);
    virtual ~ParamBase() { };

    virtual ParamBase * clone() const = 0;

    virtual void save(const boost::filesystem::path & path) const;

    virtual void read(const boost::filesystem::path & path);

    // initializes pars = {0,...}, beg = pars.begin(), end = pars.end()
    // calls initInternal, initParsSize
    virtual void init(const FixedData & fD);

    void setOffset(const unsigned int offset);

    // retrieve pars
    std::vector<double> getPar() const;
    std::vector<double> getPar(const std::vector<std::string> & name) const;

    // change pars
    // requires that newParInt has at least parsSize elements left
    // requires that beg,end point to begin() and end() of pars
    std::vector<double>::const_iterator
    putPar(std::vector<double>::const_iterator newParIt);

    // set par by name
    bool setPar(const std::string & name, const double & val);
    bool setPar(const std::vector<std::string> & name, const double & val);

    // set the probs matrix
    // adds in the corresponding term
    virtual void setFill(std::vector<double> & probs,
            const SimData & sD,
            const TrtData & tD,
            const FixedData & fD,
            const DynamicData & dD) = 0;

    // update the probs matrix for changes in simulation data
    // not an update for changes in parameters
    virtual void modFill(std::vector<double> & probs,
            const SimData & sD,
            const TrtData & tD,
            const FixedData & fD,
            const DynamicData & dD) = 0;

    virtual unsigned int size() const;

    virtual std::vector<double> partial(const int notNode,
            const int infNode,
            const SimData & sD,
            const TrtData & tD,
            const FixedData & fD,
            const DynamicData & dD) = 0;

    // second derivatives are often zero, so this is not an abstract
    // function as it will rarely need to be altered
    virtual std::vector<double> partial2(const int notNode,
            const int infNode,
            const SimData & sD,
            const TrtData & tD,
            const FixedData & fD,
            const DynamicData & dD);

    // linear scaling of the parameters
    virtual void linScale(const double & scale);
};


#endif
