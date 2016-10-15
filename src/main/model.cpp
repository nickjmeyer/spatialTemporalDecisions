#include "model.hpp"
#include <glog/logging.h>
#include "timer.hpp"

ModelBase::ModelBase(const std::string & str,
        const std::vector<ParamBase *> & newPars,
        const FixedData & fD)
    : ModelBase(str,newPars) {
    init(fD);
}


ModelBase::ModelBase(const std::string & str,
        const std::vector<ParamBase *> & newPars) {
    name = str;

    setType(INVALID);

    setFixSample(0);

    setEdgeToEdge(false);

    set = 0;
    ready = 0;
    pars = newPars;

    numPars = -1;
}


void ModelBase::init(const FixedData & fD) {
    std::for_each(pars.begin(),pars.end(),
            [&fD](ParamBase * p){
                p->init(fD);
            });

    numPars = 0;
    std::for_each(pars.begin(),pars.end(),
            [this](ParamBase * p){
                p->setOffset(this->numPars);
                this->numPars += p->size();
            });

    std::for_each(pars.begin(),pars.end(),
            [this](ParamBase * p){
                p->setTotNumPars(this->numPars);
            });
}


ModelBase::~ModelBase(){
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        delete pars.at(i);
    }
}


void ModelBase::read(){
    boost::filesystem::path paramDir = njm::sett.srcExt("");
    read_from(paramDir);
}


void ModelBase::read_from(const boost::filesystem::path path){
    boost::filesystem::path paramDir = path / ("Param"+name);
    CHECK(boost::filesystem::exists(paramDir));
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        pars[i]->read(paramDir);
    }
}


void ModelBase::save_to(const boost::filesystem::path path) const{
    boost::filesystem::path paramDir = path / ("Param"+name);
    if(!boost::filesystem::exists(paramDir)) {
        boost::system::error_code error;
        boost::filesystem::create_directories(paramDir,error);
        CHECK_EQ(error,boost::system::errc::success)
            << "failed to create directory " << paramDir;
    }

    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        pars[i]->save(paramDir);
    }
}


void ModelBase::save() const{
    boost::filesystem::path paramDir = njm::sett.srcExt("");
    save_to(paramDir);
}


void ModelBase::linScale(const double & scale){
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i)
        pars.at(i)->linScale(scale);
    set = 0;
    ready = 0;
}


void ModelBase::setType(const Estimation & est){
    fitType = est;
}


void ModelBase::setFixSample(const int & fix){
    fixSample = fix;
}

int ModelBase::getFixSample() {
    return fixSample;
}


Estimation ModelBase::getType() const{
    return fitType;
}


Estimation & ModelBase::getType() {
    return fitType;
}

void ModelBase::setEdgeToEdge(const bool edgeToEdge){
    this->set = 0;
    this->ready = 0;
    this->edgeToEdge = edgeToEdge;
}

bool ModelBase::getEdgeToEdge() const {
    return this->edgeToEdge;
}


void ModelBase::infProbs(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("infProbs");
    if(ready == 1){
        expitInfProbs.resize(sD.numNotInfec);
        int i,j,k;
        double prob;
        for(i = 0, k = 0; i < sD.numNotInfec; ++i){
            prob=1.0;
            for(j = 0; j < sD.numInfected; ++j,++k)
                prob *= quick[k];
            expitInfProbs[i] = 1.0-prob;
        }
    }
    else if(ready == 0){
        expitInfProbs.resize(sD.numNotInfec);
        int i,j;
        double prob;
        for(i = 0; i < sD.numNotInfec; ++i){
            const int iNode = sD.notInfec[i];
            prob=1.0;
            for(j = 0; j < sD.numInfected; ++j){
                const int jNode = sD.infected[j];
                if((this->getEdgeToEdge() && fD.network.at(iNode*fD.numNodes + jNode))
                        || !this->getEdgeToEdge()) {
                    prob *= 1.0 / (1.0 + std::exp(probs[iNode*fD.numNodes + jNode]));
                }
            }
            expitInfProbs[i] = 1.0-prob;
        }
    }
    // njm::timer.stop("infProbs");
}


std::vector<double> ModelBase::infProbs(){
    return expitInfProbs;
}


void ModelBase::revProbs(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("revProbs");
    if(ready == 1){
        expitRevProbs.resize(sD.numInfected);
        int i,j,k;
        double prob;
        for(i = 0,k = 0; i < sD.numInfected; ++i){
            prob=1.0;
            for(j = 0; j < sD.numNotInfec; ++j,++k)
                prob *= quick[k];
            expitRevProbs[i] = 1.0-prob;
        }
    }
    else if(ready == 0){
        expitRevProbs.resize(sD.numInfected);
        int i,j;
        double prob;
        for(i = 0; i < sD.numInfected; ++i){
            const int iNode = sD.infected[i];
            prob=1.0;
            for(j = 0; j < sD.numNotInfec; ++j) {
                const int jNode = sD.notInfec[j];
                if((this->getEdgeToEdge() && fD.network.at(iNode*fD.numNodes + jNode))
                        || !this->getEdgeToEdge()){
                    prob *= 1.0 / (1.0 + std::exp(probs[iNode*fD.numNodes + jNode]));
                }
            }
            expitRevProbs[i] = 1.0-prob;
        }
    }
    else{
        throw(1);
    }
    // njm::timer.stop("revProbs");
}


std::vector<double> ModelBase::revProbs(){
    return expitRevProbs;
}



void ModelBase::setFill(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD,
        const bool setPcPartial){
    // njm::timer.start("setFill");
    int i,parsSize = pars.size();
    probs = std::vector<double>(fD.numNodes*fD.numNodes,0.0);

    if (setPcPartial) {
        this->pcPartial.resize(fD.numNodes*fD.numNodes*numPars);
        std::fill(this->pcPartial.begin(),this->pcPartial.end(),0.);

        for(i = 0; i < parsSize; ++i){
            pars[i]->setFill(probs,this->pcPartial,sD,tD,fD,dD);
        }
    } else {
        for(i = 0; i < parsSize; ++i){
            pars[i]->setFill(probs,sD,tD,fD,dD);
        }
    }

    set = 1;
    ready = 0;
    // njm::timer.stop("setFill");
}


void ModelBase::modFill(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD,
        const bool setPcPartial){
    // njm::timer.start("modFill");
    if(set == 1){
        int i,parsSize = pars.size();

        if (setPcPartial) {
            for(i = 0; i < parsSize; ++i){
                pars[i]->modFill(probs,pcPartial,sD,tD,fD,dD);
            }
        } else {
            for(i = 0; i < parsSize; ++i){
                pars[i]->modFill(probs,sD,tD,fD,dD);
            }
        }
    }
    else if(set == 0){
        setFill(sD,tD,fD,dD,setPcPartial);
    }
    ready = 0;
    // njm::timer.stop("modFill");
}


void ModelBase::setQuick(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("setQuick");
    int i,j,k,pK;
    quick.resize(sD.numNotInfec * sD.numInfected);
    for(i = 0,k = 0; i < sD.numNotInfec; ++i){
        const int iNode = sD.notInfec[i];
        pK = iNode*fD.numNodes;
        for(j = 0; j < sD.numInfected; ++j,++k){
            const int jNode = sD.infected[j];
            if((this->getEdgeToEdge() && fD.network.at(iNode*fD.numNodes + jNode))
                    || !this->getEdgeToEdge()){
                quick[k] = 1.0/(1.0 + std::exp(probs[pK + jNode]));
            } else {
                quick[k] = 1.0;
            }
        }
    }
    ready = 1;
    // njm::timer.stop("setQuick");
}


std::vector<double> & ModelBase::getQuick() {
    return quick;
}



double ModelBase::oneOnOne(const int notNode,
        const int infNode,
        const int numNodes) const {
    return probs[notNode * numNodes + infNode];
}



std::vector<double> ModelBase::getPar() const{
    std::vector<double> all,add;
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        add = pars[i]->getPar();
        all.insert(all.end(),add.begin(),add.end());
    }
    return all;
}


std::vector<double>
ModelBase::getPar(const std::vector<std::string> & name) const{
    std::vector<double> res,add;
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        add = pars[i]->getPar(name);
        res.insert(res.end(),add.begin(),add.end());
    }
    return res;
}


std::vector<double>::const_iterator
ModelBase::putPar(std::vector<double>::const_iterator it){
    set = 0;
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        it = pars[i]->putPar(it);
    }
    return it;
}


void ModelBase::setPar(const std::string & name,
        const double & val){
    int i,parsSize = pars.size();
    bool found = false;
    for(i = 0; i < parsSize; ++i){
        found |= pars[i]->setPar(name,val);
    }

    if(!found){
        std::cout << "name " + name + " not found in setPar()" << std::endl;
        throw(1);
    }
}


void ModelBase::setPar(const std::vector<std::string> & name,
        const double & val){
    int i,parsSize = pars.size();
    bool found = false;
    for(i = 0; i < parsSize; ++i)
        found |= pars[i]->setPar(name,val);

    if(!found){
        std::cout << "name not found in setPar()" << std::endl;
        throw(1);
    }
}


std::vector<double> ModelBase::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("partial");
    std::vector<double> p,pi;
    p.reserve(numPars);
    int i,parsSize = pars.size();
    for(i = 0; i < parsSize; ++i){
        pi = pars[i]->partial(notNode,infNode,sD,tD,fD,dD);
        p.insert(p.end(),pi.begin(),pi.end());
    }

    // njm::timer.stop("partial");
    return p;
}


std::vector<double> ModelBase::partial2(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("partial2");
    std::vector<double> p,pi;
    std::vector<int> parsLen;
    int i,parsSize = pars.size(),totLen = 0;
    for(i = 0; i < parsSize; ++i){
        parsLen.push_back(totLen);
        totLen += pars[i]->size();
    }
    parsLen.push_back(totLen);


    p.resize(totLen*totLen);
    std::fill(p.begin(),p.end(),0);

    int j,k,n,N,D;
    for(i = 0; i < parsSize; ++i){
        pi = pars[i]->partial2(notNode,infNode,sD,tD,fD,dD);

        n = parsLen[i];
        N = parsLen[i+1];
        D = N - n;
        for(j = n; j < N; ++j){
            for(k = j; k < N; ++k){
                p[j*totLen + k] = pi[(j-n)*D + (k-n)];
                p[k*totLen + j] = pi[(k-n)*D + (j-n)];
            }
        }

    }

    // njm::timer.stop("partial2");
    return p;
}


void ModelBase::setFisher(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("setFisher");
    std::vector<double> hess = logllHess(sD,tD,fD,dD);
    arma::mat I(hess.data(),numPars,numPars);

    if(sD.time <= fD.trtStart){
        // if no treatments, then can't estimate asymptotic distribution
        // set asymptotic variance as 1.0
        // assign -1.0 since I gets scaled by -1.0 later
        I.row(numPars-2) = arma::zeros<arma::rowvec>(numPars);
        I.row(numPars-1) = arma::zeros<arma::rowvec>(numPars);
        I.col(numPars-2) = arma::zeros<arma::colvec>(numPars);
        I.col(numPars-1) = arma::zeros<arma::colvec>(numPars);

        I(numPars-2,numPars-2) = -1.0;
        I(numPars-1,numPars-1) = -1.0;
    }

    // information matrix is - E[ d^2/(d\theta)^2 log L(\theta) ]
    I *= -1.0;

    arma::mat eigvec;
    arma::colvec eigval;
    arma::eig_sym(eigval,eigvec,I);

    std::vector<double> currPar = getPar();
    meanHit = arma::colvec(currPar.data(),numPars);


    // invert the non-zero eigen values
    unsigned int pi;
    for(pi = 0; pi < numPars; ++pi){
        if(eigval(pi) < 1e-10)
            // stable non-negative definite matrix
            eigval(pi) = 0.0;
        else{
            eigval(pi) = 1.0/std::sqrt(eigval(pi));
        }
    }
    varHit = eigvec * arma::diagmat(eigval);
    // njm::timer.stop("setFisher");
}


bool ModelBase::sample(const bool force){
    if(fitType == MLES && (!fixSample || force)){
        arma::colvec rand(numPars);
        unsigned int i;
        for(i = 0; i < numPars; ++i){
            rand(i) = njm::rnorm01();
        }

        std::vector<double> sample =
            arma::conv_to<std::vector<double> >::from(meanHit + varHit * rand);

        putPar(sample.begin());

        return true;
    }
    else{
        return false;
    }
}


void ModelBase::revert(){
    if(fitType == MLES){

        std::vector<double> sample =
            arma::conv_to<std::vector<double> >::from(meanHit);

        putPar(sample.begin());
    }
    else{
        std::cout << "ModelBase::revert(): not implemented for fytType of "
                  << fitType << std::endl;
        throw(1);
    }
}



void ModelBase::fit(const SimData & sD, const TrtData & tD,
        const FixedData & fD, const DynamicData & dD,
        const bool init){
    if(init){
        // used penalized version of current par;
        std::vector<double> currPar = getPar();
        for (int i = 0; i < currPar.size(); ++i) {
            currPar.at(i) *= 0.25;
        }
        fit(currPar,sD,tD,fD,dD);
    }
    else{
        std::vector<double> startingVals(numPars,0.2);
        startingVals.at(0) = -3.0;
        fit(startingVals,sD,tD,fD,dD);
    }
}


void ModelBase::fit(const SimData & sD, const TrtData & tD,
        const FixedData & fD, const DynamicData & dD){

    std::vector<double> startingVals(numPars,0.2);
    startingVals.at(0) = -3.0;
    fit(startingVals,sD,tD,fD,dD);
}


void ModelBase::fit(const std::vector<double> & startingVals,
        const SimData & sD, const TrtData & tD,
        const FixedData & fD, const DynamicData & dD){
    if(fitType == MLE || fitType == MLES){

        // estimate the MLE
        estimateMle(startingVals,sD,tD,fD,dD);


        if(sD.time <= fD.trtStart){
            // if before trt starts, set trt values to prior
            std::vector<std::string> names = {"trtPre","trtAct"};
            setPar(names,fD.priorTrtMean);
        }


        if(fitType == MLES){
            // for sampling
            setFisher(sD,tD,fD,dD);
        }

        // for probabilities, must do after fisher
        setFill(sD,tD,fD,dD);

    }
    else if(fitType == MCMC){
        std::cout << "Error: ModelBase::fit(): MCMC not setup"
                  << std::endl;
        throw(1);
    }
    else{
        std::cout << "Not a valid Estimation of "
                  << fitType
                  << std::endl;
        throw(1);
    }
}



void ModelBase::estimateMle(const std::vector<double> & startingVals,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){

    if(fitType != MLE && fitType != MLES){
        std::cout << "ModelBase::estimateMLE is called for non MLE or MLES type"
                  << std::endl;
        throw(1);
    }

    ModelBaseFitObj fitObj(this,sD,tD,fD,dD);

    const gsl_multimin_fdfminimizer_type * T;
    gsl_multimin_fdfminimizer *s;

    gsl_vector * x;
    x = gsl_vector_alloc(this->numPars);
    int pi;
    for(pi = 0; pi < int(this->numPars); ++pi){
        gsl_vector_set(x,pi,startingVals.at(pi));
    }

    gsl_multimin_function_fdf my_func;
    my_func.n = this->numPars;
    my_func.f = objFn;
    my_func.df = objFnGrad;
    my_func.fdf = objFnBoth;
    my_func.params = &fitObj;

    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T,this->numPars);

    gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);

    int iter = 0;
    int status;
    const int maxIter = 100;
    do{
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if(status)
            break;

        status = gsl_multimin_test_gradient(s->gradient,0.1);

    }while(status == GSL_CONTINUE && iter < maxIter);

    std::vector<double> gradVals;
    for (pi = 0; pi < int(numPars); ++pi) {
        gradVals.push_back(gsl_vector_get(s->gradient,pi));
    }

    // CHECK_LT(iter,maxIter)
    //     << "Reached max iterations"
    //     << std::endl
    //     << "status: " << status << std::endl
    //     << "iter: " << iter << std::endl
    //     << "numInfected: " << sD.numInfected << std::endl
    //     << "numNotInfec: " << sD.numNotInfec << std::endl
    //     << "time: " << sD.time << std::endl
    //     << "gradient check: "
    //     << gsl_multimin_test_gradient(s->gradient,0.1) << std::endl
    //     << "gradient: " << njm::toString(gradVals," ","") << std::endl;
    // if(iter >= maxIter)
    //     std::cout << "exceeded iters" << std::endl;

    CHECK(status == GSL_SUCCESS || status == GSL_CONTINUE)
        << std::endl
        << "status: " << status << std::endl
        << "iter: " << iter << std::endl
        << "numInfected: " << sD.numInfected << std::endl
        << "time: " << sD.time << std::endl
        << "gradient check: "
        << gsl_multimin_test_gradient(s->gradient,0.1) << std::endl
        << "f: " << s->f << std::endl
        << "gradient: " << njm::toString(gradVals," ","") << std::endl
        << "starting: " << njm::toString(startingVals," ","") << std::endl;

    // CHECK(status == GSL_SUCCESS || status == GSL_CONTINUE ||
    //         (status == 27
    //                 && sD.numInfected == 1
    //                 && sD.time > fD.trtStart))
    //     << std::endl
    //     << "status: " << status << std::endl
    //     << "iter: " << iter << std::endl
    //     << "numInfected: " << sD.numInfected << std::endl
    //     << "time: " << sD.time << std::endl
    //     << "gradient check: "
    //     << gsl_multimin_test_gradient(s->gradient,0.1) << std::endl
    //     << "f: " << s->f << std::endl
    //     << "gradient: " << njm::toString(gradVals," ","") << std::endl
    //     << "starting: " << njm::toString(startingVals," ","") << std::endl;

// #pragma omp critical
//     {
//         if(status != GSL_SUCCESS) {
//             std::cout << "setup values" << std::endl;
//             std::cout << "time: " << sD.time << std::endl
//                       << "numInfected: " << sD.numInfected << std::endl
//                       << "numNotInfec: " << sD.numNotInfec << std::endl
//                       << "sum infected: " << std::accumulate(sD.infected.begin(),
//                               sD.infected.end(),0) << std::endl
//                       << "sum notInfec: " << std::accumulate(sD.notInfec.begin(),
//                               sD.notInfec.end(),0) << std::endl;
//             for (int i = 0; i < sD.history.size(); ++i) {
//                 std::cout << "sum history " << i << ": "
//                           << std::accumulate(sD.history.at(i).begin(),
//                                   sD.history.at(i).end(),0)
//                           << std::endl;
//             }

//             std::cout << "Gradient: "
//                       << njm::toString(gradVals," ","") << std::endl;
//             std::cout << "Status: " << status << std::endl;
//             CHECK_LT(iter,maxIter) << "Reached maximum iterations";
//             CHECK_EQ(status,GSL_SUCCESS)
//                 << std::endl
//                 << "Iterations: " << iter << std::endl
//                 << "Num infected: " << sD.numInfected << std::endl
//                 << "Infected: " << njm::toString(sD.infected," ","") << std::endl
//                 << "Num notInfec: " << sD.numNotInfec << std::endl
//                 << "Gradient: " << njm::toString(gradVals," ","") << std::endl;
//         }
//     }


    std::vector<double> mle;
    for(pi = 0; pi < int(numPars); ++pi){
        mle.push_back(gsl_vector_get(s->x,pi));
    }

    this->putPar(mle.begin());


    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
}



double ModelBase::logll(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){

    std::vector<std::vector<int> > hist = sD.history;
    hist.push_back(sD.status);
    std::vector<DataBundle> db = historyToData(hist);

    double logllVal = 0.0;

    int t,nN;
    // loop over time points
    setFill(sD,tD,fD,dD);
    for(t = 0; t < sD.time; ++t){
        const SimData & sDi = std::get<0>(db[t]);
        const TrtData & tDi = std::get<1>(db[t]);
        const DynamicData & dDi = std::get<2>(db[t]);

        modFill(sDi,tDi,fD,dDi);
        infProbs(sDi,tDi,fD,dDi);

        if(int(expitInfProbs.size()) != sDi.numNotInfec){
            std::cout << "ModelBase::logll(): length of expitInfProbs is not same as"
                      << " number of uninfected nodes at time t"
                      << std::endl;
            throw(1);
        }

        // njm::timer.start("logll_computation");
        // loop over uninfected nodes at time t
        for(nN = 0; nN < sDi.numNotInfec; ++nN){
            double prob = expitInfProbs.at(nN);
            int next = (hist[t+1][sDi.notInfec[nN]] < 2) ? 0 : 1;
            if(next == 1){
                if(prob < 1e-44)
                    logllVal += -100.0;
                else
                    logllVal += std::log(prob);
            }
            else{
                if((1.0-prob) < 1e-44)
                    logllVal += -100.0;
                else
                    logllVal += std::log(1.0 - prob);
            }
        }
        // njm::timer.stop("logll_computation");
    }
    return logllVal;
}


std::vector<double> ModelBase::logllGrad(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    std::vector<std::vector<int> > hist = sD.history;
    hist.push_back(sD.status);
    std::vector<DataBundle> db = historyToData(hist);

    std::vector<double> logllGradVal(numPars,0.0);
    std::fill(logllGradVal.begin(),logllGradVal.end(),0.0);

    setFill(sD,tD,fD,dD,true);
    int t,nN,iN,pi;
    // loop over time points
    for(t = 0; t < sD.time; ++t){
        const SimData & sDi = std::get<0>(db[t]);
        const TrtData & tDi = std::get<1>(db[t]);
        const DynamicData & dDi = std::get<2>(db[t]);

        modFill(sDi,tDi,fD,dDi,true);
        setQuick(sDi,tDi,fD,dDi);
        infProbs(sDi,tDi,fD,dDi);

        if(int(expitInfProbs.size()) != sDi.numNotInfec){
            std::cout << "ModelBase::logll(): length of "
                      << "expitInfProbs is not same as"
                      << " number of uninfected nodes at time t"
                      << std::endl;
            throw(1);
        }

        // njm::timer.start("logllGrad_computation");
        // loop over uninfected nodes at time t
        for(nN = 0; nN < sDi.numNotInfec; ++nN){
            const int nNode = sDi.notInfec[nN];
            double prob = expitInfProbs.at(nN);
            int next = (hist[t+1][nNode] < 2) ? 0 : 1;

            if(prob > 0.0){
                double beg = double(next)/prob - 1.0;
                for(iN = 0; iN < sDi.numInfected; ++iN){
                    const int iNode = sDi.infected[iN];
                    if((this->getEdgeToEdge()
                                    && fD.network.at(nNode*fD.numNodes + iNode))
                            || !this->getEdgeToEdge()){
                        // Quick stores probabilty of not infecting.
                        // Need to take (1.0 - quick) to get
                        // probablity of infecting.

                        double nNInfByiN = 1.0 - quick[nN*sDi.numInfected + iN];
                        const int partialIndex = nNode*fD.numNodes*numPars +
                            iNode*numPars;
                        for(pi=0; pi < int(numPars); ++pi){
                            logllGradVal.at(pi) += beg * nNInfByiN *
                                pcPartial.at(partialIndex + pi);
                        }
                    }
                }
            }
        }
        // njm::timer.stop("logllGrad_computation");
    }
    return logllGradVal;
}



std::pair<double, std::vector<double> > ModelBase::logllBoth(
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){

    std::vector<std::vector<int> > hist = sD.history;
    hist.push_back(sD.status);
    std::vector<DataBundle> db = historyToData(hist);

    std::vector<double> logllGradVal(numPars,0.0);
    std::fill(logllGradVal.begin(),logllGradVal.end(),0.0);

    double logllVal = 0.0;

    setFill(sD,tD,fD,dD,true);
    int t,nN,iN,pi;
    // loop over time points
    for(t = 0; t < sD.time; ++t){
        const SimData & sDi = std::get<0>(db[t]);
        const TrtData & tDi = std::get<1>(db[t]);
        const DynamicData & dDi = std::get<2>(db[t]);

        modFill(sDi,tDi,fD,dDi,true);
        setQuick(sDi,tDi,fD,dDi);
        infProbs(sDi,tDi,fD,dDi);

        if(int(expitInfProbs.size()) != sDi.numNotInfec){
            std::cout << "ModelBase::logll(): length of expitInfProbs "
                      << "is not same as number of uninfected nodes at time t"
                      << std::endl;
            throw(1);
        }


        // loop over uninfected nodes at time t
        for(nN = 0; nN < sDi.numNotInfec; ++nN){
            // njm::timer.start("logll_computation");
            const int nNode = sDi.notInfec[nN];
            const double prob = expitInfProbs.at(nN);
            const int next = (hist[t+1][nNode] < 2) ? 0 : 1;

            // log likelihood
            if(next == 1){
                if(prob < 1e-44)
                    logllVal += -100.0;
                else
                    logllVal += std::log(prob);
            }
            else{
                if((1.0-prob) < 1e-44)
                    logllVal += -100.0;
                else
                    logllVal += std::log(1.0 - prob);
            }
            // njm::timer.stop("logll_computation");

            // njm::timer.start("logllGrad_computation");
            // log likelihood gradient
            if(prob > 0.0){
                const double beg = double(next)/prob - 1.0;
                for(iN = 0; iN < sDi.numInfected; ++iN){
                    const int iNode = sDi.infected[iN];
                    if((this->getEdgeToEdge() &&
                                    fD.network.at(nNode*fD.numNodes + iNode))
                            || !this->getEdgeToEdge()){
                        // Quick stores probabilty of not infecting.
                        // Need to take (1.0 - quick) to get
                        // probablity of infecting.

                        const double nNInfByiN =
                            1.0 - quick[nN*sDi.numInfected + iN];
                        const int partialIndex = nNode*fD.numNodes*numPars
                            + iNode*numPars;
                        for(pi=0; pi < int(numPars); ++pi){
                            logllGradVal.at(pi) += beg * nNInfByiN *
                                pcPartial.at(partialIndex + pi);
                        }
                    }
                }
            }
            // njm::timer.stop("logllGrad_computation");
        }
    }
    return std::pair<double,std::vector<double> >(logllVal,logllGradVal);
}



std::vector<double> ModelBase::logllHess(const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    std::vector<double> hess(numPars * numPars);
    std::fill(hess.begin(),hess.end(),0);

    std::vector<std::vector<int> > hist = sD.history;
    hist.push_back(sD.status);
    std::vector<DataBundle> db = historyToData(hist);

    int t,iN,nN;
    unsigned int pi,pj;
    double addVal = 0.0;
    for(t = 0; t < sD.time; ++t){
        SimData sDi= std::get<0>(db[t]);
        TrtData tDi= std::get<1>(db[t]);
        DynamicData dDi = std::get<2>(db[t]);

        setFill(sDi,tDi,fD,dDi);
        setQuick(sDi,tDi,fD,dDi);
        infProbs(sDi,tDi,fD,dDi);

        for(nN = 0; nN < sDi.numNotInfec; ++nN){
            const int nNode = sDi.notInfec[nN];
            std::vector<double> dbl(numPars*numPars,0);
            std::vector<double> sqr(numPars,0);

            for(iN = 0; iN < sDi.numInfected; ++iN){
                const int iNode = sDi.infected[iN];
                if((this->getEdgeToEdge() && fD.network.at(nNode*fD.numNodes + iNode))
                        || !this->getEdgeToEdge()){
                    int ind = nN*sDi.numInfected + iN;
                    std::vector<double> p = partial(nNode,iNode,
                            sDi,tDi,fD,dDi);
                    std::vector<double> p2 = partial2(nNode,iNode,
                            sDi,tDi,fD,dDi);

                    // double quickInd = std::max(1e-10,quick[ind]);
                    double quickInd = 1.0 - quick[ind];
                    for(pi = 0; pi < numPars; ++pi){
                        sqr[pi] += quickInd*p[pi];

                        int piInd = pi*numPars;
                        for(pj = pi; pj < numPars; ++pj){
                            addVal = quickInd*p2[piInd + pj];
                            addVal += quickInd*(1.0-quickInd)*p[pj]*p[pi];

                            dbl[piInd + pj] += addVal;
                            if(pj != pi){
                                dbl[pj*numPars + pi] += addVal;
                            }
                        }
                    }
                }
            }

            int next = (hist[t+1][sDi.notInfec[nN]] < 2 ? 0 : 1);

            for(pi = 0; pi < numPars; ++pi){
                for(pj = pi; pj < numPars; ++pj){
                    // double prob = std::max(1e-10,expitInfProbs[nN]);
                    double prob = std::max(1e-10,expitInfProbs[nN]);

                    if(prob > 0.0){
                        addVal = (double(next)/prob - 1.0)*dbl[pi*numPars + pj];
                        hess[pi*numPars + pj] += addVal;

                        if(pj != pi){
                            hess[pj*numPars + pi] += addVal;
                        }
                    }

                    if((prob*prob) > 0.0){
                        addVal = (1-prob)*sqr[pi]*sqr[pj]/(prob*prob);

                        if(next == 1){
                            hess[pi*numPars + pj] -= addVal;

                            if(pj != pi){
                                hess[pj*numPars + pi] -= addVal;
                            }
                        }
                    }
                }
            }
        }
    }

    return hess;
}



ModelBaseFitObj::ModelBaseFitObj(ModelBase * const mb,
        const SimData sD,
        const TrtData tD,
        const FixedData fD,
        const DynamicData dD){
    this->mb = mb;
    this->sD = sD;
    this->tD = tD;
    this->fD = fD;
    this->dD = dD;
}


double objFn(const gsl_vector * x, void * params){
    // njm::timer.start("objFn");
    ModelBaseFitObj * fitObj = static_cast<ModelBaseFitObj*>(params);
    std::vector<double> par;
    int pi;
    for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
        par.push_back(gsl_vector_get(x,pi));
    }

    fitObj->mb->putPar(par.begin());

    // return negative since GSL minimizes the function
    double ll = fitObj->mb->logll(fitObj->sD,fitObj->tD,fitObj->fD,fitObj->dD);
    CHECK(std::isfinite(ll)) << "Likelihood value is not finite";
    // njm::timer.stop("objFn");
    return - ll;
}

void objFnGrad(const gsl_vector * x, void * params, gsl_vector * g){
    // njm::timer.start("objFnGrad");
    ModelBaseFitObj * fitObj = static_cast<ModelBaseFitObj*>(params);
    std::vector<double> par;
    int pi;
    for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
        par.push_back(gsl_vector_get(x,pi));
    }

    fitObj->mb->putPar(par.begin());

    std::vector<double> llGrad = fitObj->mb->logllGrad(fitObj->sD,fitObj->tD,
            fitObj->fD,fitObj->dD);
    for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
        // assign the negative of the gradient value
        // GSL minimizes the function, need to adjust the gradient too
        CHECK(std::isfinite(llGrad.at(pi)))
            << "Likelihood gradient value is not finite for parameter index "
            << pi << " for model " << fitObj->mb->name;
        gsl_vector_set(g,pi,-llGrad.at(pi));
    }
    // njm::timer.stop("objFnGrad");
}


void objFnBoth(const gsl_vector * x, void * params, double * f, gsl_vector * g){
    // njm::timer.start("objFnBoth");
    ModelBaseFitObj * fitObj = static_cast<ModelBaseFitObj*>(params);
    std::vector<double> par;
    int pi;
    for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
        par.push_back(gsl_vector_get(x,pi));
    }

    fitObj->mb->putPar(par.begin());

    std::pair<double,std::vector<double> > both =
        fitObj->mb->logllBoth(fitObj->sD,fitObj->tD,
                fitObj->fD,fitObj->dD);

    // log ll
    CHECK(std::isfinite(both.first));
    *f = -both.first;

    // log ll grad
    for(pi = 0; pi < int(fitObj->mb->numPars); ++pi){
        // assign the negative of the gradient value
        // GSL minimizes the function, need to adjust the gradient too
        CHECK(std::isfinite(both.second.at(pi)))
            << "Likelihood gradient value is not finite for parameter index "
            << pi << " for model " << fitObj->mb->name;
        gsl_vector_set(g,pi,-both.second.at(pi));
    }
    // njm::timer.stop("objFnBoth");
}
