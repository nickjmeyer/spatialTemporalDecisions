#include "paramTrt.hpp"
#include "timer.hpp"



unsigned int ParamTrt::initParsSize(const FixedData & fD){
    return 2U;
}


std::vector<std::string> ParamTrt::initNames(){
    return {"trtAct","trtPre"};
}


void ParamTrt::initInternal(const FixedData & fD){
    numNodes = fD.numNodes;

    a = std::vector<int>(numNodes,0);
    p = std::vector<int>(numNodes,0);
}


void ParamTrt::updateBefore(){
}


void ParamTrt::updateAfter(){
}


void ParamTrt::setFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    double trtAct = pars.at(0);
    double trtPre = pars.at(1);

    int i,j,k;
    a = tD.a;
    p = tD.p;
    for(i = 0; i < numNodes; ++i){
        if(p.at(i) == 1){ // add if pre
            k = i*numNodes;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j) -= trtPre;
            }
        }
        if(a.at(i) == 1){ // add if act
            k = i;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j*numNodes) -= trtAct;
            }
        }
    }
}


void ParamTrt::setFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    double trtAct = pars.at(0);
    double trtPre = pars.at(1);

    int i,j,k;
    a = tD.a;
    p = tD.p;
    for(i = 0; i < numNodes; ++i){
        if(p.at(i) == 1){ // add if pre
            k = i*numNodes;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j) -= trtPre;
                pcPartial.at(k*totNumPars + j*totNumPars) = -trtPre;
            }
        }
        if(a.at(i) == 1){ // add if act
            k = i;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j*numNodes) -= trtAct;
                pcPartial.at(k*numNodes*totNumPars + j*totNumPars) = -trtAct;
            }
        }
    }
}


void ParamTrt::modFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    const double trtAct = pars.at(0);
    const double trtPre = pars.at(1);
    int i,j,k;
    for(i = 0; i < numNodes; ++i){
        if(tD.p.at(i) == 1 && p.at(i) == 0){ // add pre
            k = i*numNodes;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j) -= trtPre;
            }
        }
        else if(tD.p.at(i) == 0 && p.at(i) == 1){ // remove pre
            k = i*numNodes;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j) += trtPre;
            }
        }
    }

    for(i = 0; i < numNodes; ++i){
        if(tD.a.at(i) == 1 && a.at(i) == 0){ // add act
            k = i;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j*numNodes) -= trtAct;
            }
        }
        if(tD.a.at(i) == 0 && a.at(i) == 1){ // remove act
            k = i;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j*numNodes) += trtAct;
            }
        }
    }

    a = tD.a;
    p = tD.p;
}


void ParamTrt::modFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    const double trtAct = pars.at(0);
    const double trtPre = pars.at(1);
    int i,j,k;
    for(i = 0; i < numNodes; ++i){
        if(tD.p.at(i) == 1 && p.at(i) == 0){ // add pre
            k = i*numNodes;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j) -= trtPre;
                // pre is second, hence offset + 1
                pcPartial.at(k*totNumPars + j*totNumPars + offset + 1) -= 1;
            }
        }
        else if(tD.p.at(i) == 0 && p.at(i) == 1){ // remove pre
            k = i*numNodes;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j) += trtPre;
                // pre is second, hence offset + 1
                pcPartial.at(k*totNumPars + j*totNumPars + offset + 1) += 1;
            }
        }
    }

    for(i = 0; i < numNodes; ++i){
        if(tD.a.at(i) == 1 && a.at(i) == 0){ // add act
            k = i;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j*numNodes) -= trtAct;
                pcPartial.at(k*totNumPars + j*numNodes*totNumPars +
                        offset) -= 1; // act is first, hence offset + 0
            }
        }
        if(tD.a.at(i) == 0 && a.at(i) == 1){ // remove act
            k = i;
            for(j = 0; j < numNodes; ++j){
                probs.at(k + j*numNodes) += trtAct;
                pcPartial.at(k*totNumPars + j*numNodes*totNumPars +
                        offset) += 1; // act is first, hence offset + 0
            }
        }
    }

    a = tD.a;
    p = tD.p;
}


std::vector<double> ParamTrt::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    njm::timer.start("partial_trt");
    std::vector<double> p;
    p.push_back(-tD.a.at(infNode));
    p.push_back(-tD.p.at(notNode));
    njm::timer.stop("partial_trt");
    return p;
}
