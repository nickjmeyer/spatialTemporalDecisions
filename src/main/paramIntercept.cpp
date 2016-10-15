#include "paramIntercept.hpp"


unsigned int ParamIntercept::initParsSize(const FixedData & fD){
    return 1U; // unsigned int literal
}

std::vector<std::string> ParamIntercept::initNames(){
    return {"intcp"};
}


void ParamIntercept::initInternal(const FixedData & fD){
}


void ParamIntercept::updateBefore(){
}


void ParamIntercept::updateAfter(){
}


void ParamIntercept::setFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    const double val = *beg;
    std::for_each(probs.begin(),probs.end(),
            [&val](double & x){
                x += val;
            });
}


void ParamIntercept::setFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    const double val = *beg;
    std::for_each(probs.begin(),probs.end(),
            [&val](double & x){
                x += val;
            });

    const int I = fD.numNodes*fD.numNodes;
    for (int i = 0; i < I; i++) {
        pcPartial.at(i*totNumPars + offset) = 1;
    }
}


void ParamIntercept::modFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // do nothing
}


void ParamIntercept::modFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // do nothing
}


std::vector<double> ParamIntercept::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    return std::vector<double>(parsSize,1);
}
