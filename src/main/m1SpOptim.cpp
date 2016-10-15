#include "m1SpOptim.hpp"


M1SpOptimTunePar::M1SpOptimTunePar(){
    mcReps = 10;

    C = 10.0;

    t = 1.0;

    ell = 1.25;

    muMin = 0.1;

    A = 30;
    B = 1;

    tune = 0;

    fixSample = 0;
}


std::vector<double> M1SpOptimTunePar::getPar() const{
    std::vector<double> par = {A,B};
    return par;
}


void M1SpOptimTunePar::putPar(const std::vector<double> & par){
    std::vector<double>::const_iterator it;
    it = par.begin();
    A = *it++;
    B = *it++;
}


template <class S, class A, class M>
M1SpOptim<S,A,M>::M1SpOptim(){
    name = "M1Sp";
}


template <class S, class A, class M>
void M1SpOptim<S,A,M>::reset(){
    tp.A = 30;
    tp.B = 1;
}


template <class S, class A, class M>
void M1SpOptim<S,A,M>
::optim(const S & system,
        A & agent){

    System<M,M> s(system.sD,system.tD,system.fD,system.dD,
            system.modelEst,system.modelEst);

    if(tp.tune != 0 && system.sD.time == (system.fD.trtStart + 1))
        tune(s,agent);

    if(tp.fixSample != 0){
        s.modelGen_r.setFixSample(1);
        s.modelEst_r.setFixSample(1);

        s.revert();
    }


    PlainRunner<System<M,M>,A> runner;

    std::vector<double> par=agent.tp.getPar();
    int i,converged=0,numPar = par.size();
    std::vector<double> h(numPar,0.0);
    std::vector<double> parPH(numPar,0.0),parMH(numPar,0.0);


    double valP, valM;

    int iter=1;

    double mu = tp.A/std::pow(tp.B+iter,tp.ell);
    double cm = tp.C / std::pow(iter,tp.t);

    while(!converged){

        for(i=0; i<numPar; i++){
            h.at(i) = (njm::rber(0.5) == 1 ? 1.0 : -1.0) * cm;
            parPH.at(i) = par.at(i) + h.at(i);
            parMH.at(i) = par.at(i) - h.at(i);
        }

        if(tp.fixSample != 0){
            s.modelGen_r.sample(true);
            s.modelEst_r = s.modelGen_r;

            s.revert();
        }


        njm::timer.start("optim run");
        agent.tp.putPar(parPH);
        valP = runner.run(s,agent,tp.mcReps,s.fD.finalT).sMean();

        agent.tp.putPar(parMH);
        valM = runner.run(s,agent,tp.mcReps,s.fD.finalT).sMean();
        njm::timer.stop("optim run");


        for(i=0; i<numPar; i++)
            par.at(i) = par.at(i) - mu*(valP - valM)/(2.0*h.at(i));


        // if(omp_get_thread_num() == 0)
        //   std::cout << "iter: " + njm::toString(iter,"",4,0) +
        // 	" || " + njm::toString(valP,"",6,4) + " - " +
        // 	njm::toString(valM,"",6,4) + " -> " +
        // 	njm::toString(mu,"",6,4) + " , " + njm::toString(cm,"",6,4) +
        // 	" || " + njm::toString(par,", ","") << "\r" << std::flush;


        ++iter;

        mu = tp.A/std::pow(tp.B + iter,tp.ell);
        cm = tp.C/std::pow(iter,tp.t);


        if(mu < tp.muMin){
            converged = 1;
        }

    }

    agent.tp.putPar(par); // assign optimized par to the agent
}



template <class S, class A, class M>
void M1SpOptim<S,A,M>
::tune(const System<M,M> & system,
        A agent){

    std::cout << "thread "
              << omp_get_thread_num()
              << " is tuning!!!!!!!" << std::endl;
    System<M,M> s(system.sD_r,system.tD_r,system.fD,system.dD_r,
            system.modelEst,system.modelEst);
    s.modelEst.fitType = MLE;
    s.fD.finalT = s.sD.time + 2*s.fD.period;

    M1SpOptim<System<M,M>,A,M> o;
    o.tp.tune = 0;

    TuneRunner<System<M,M>,A,
               M1SpOptim<System<M,M>,A,M> > r;

    std::vector<double> scale;
    scale.push_back(0.5);
    scale.push_back(1.0);
    scale.push_back(2.0);


    int i,j;
    std::vector<std::pair<double,double> > abVals;
    for(i = 0; i < (int)scale.size(); ++i)
        for(j = 0; j < (int)scale.size(); ++j)
            abVals.push_back(std::pair<double,double>(30*scale.at(i),
                            1*scale.at(j)));

    int numAbVals=abVals.size();
    double val,minVal=1.0,bestA=10,bestB=100;
    for(i=0; i<numAbVals; i++){
        o.tp.A=abVals.at(i).first;
        o.tp.B=abVals.at(i).second;

        val = r.run(s,agent,o,50,s.fD.finalT).sMean();
        if(val < minVal){
            bestA = o.tp.A;
            bestB = o.tp.B;

            minVal = val;
        }
    }

    tp.A = bestA;
    tp.B = bestB;

}


template class M1SpOptim<System<ModelGDist,
                                ModelGDist>,
                         RankAgent<WnsFeatures3<ModelGDist>,
                                   ModelGDist>,
                         ModelGDist>;


template class M1SpOptim<System<ModelGravityGDist,
                                ModelGravityGDist>,
                         RankAgent<WnsFeatures3<ModelGravityGDist>,
                                   ModelGravityGDist>,
                         ModelGravityGDist>;

template class M1SpOptim<System<Model2GravityGDist,
                                Model2GravityGDist>,
                         RankAgent<WnsFeatures3<Model2GravityGDist>,
                                   Model2GravityGDist>,
                         Model2GravityGDist>;

template class M1SpOptim<System<Model2EdgeToEdge,
                                ModelIntercept>,
                         RankAgent<WnsFeatures3<ModelIntercept>,
                                   ModelIntercept>,
                         ModelIntercept>;

template class M1SpOptim<System<Model2EdgeToEdge,
                                Model2EdgeToEdge>,
                         RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                                   Model2EdgeToEdge>,
                         Model2EdgeToEdge>;

template class M1SpOptim<System<Model2GravityEDist,
                                Model2GravityEDist>,
                         RankAgent<WnsFeatures3<Model2GravityEDist>,
                                   Model2GravityEDist>,
                         Model2GravityEDist>;

template class M1SpOptim<System<Model2GravityEDist,
                                ModelEDist>,
                         RankAgent<WnsFeatures3<ModelEDist>,
                                   ModelEDist>,
                         ModelEDist>;

template class M1SpOptim<System<Model2GravityGDist,
                                ModelGDist>,
                         RankAgent<WnsFeatures3<ModelGDist>,
                                   ModelGDist>,
                         ModelGDist>;

template class M1SpOptim<System<Model2GravityGDist,
                                Model2GravityGDist>,
                         RankAgent<ToyFeatures5<Model2GravityGDist>,
                                   Model2GravityGDist>,
                         Model2GravityGDist>;

template class M1SpOptim<System<Model2GravityEDist,
                                Model2GravityEDist>,
                         RankAgent<ToyFeatures5<Model2GravityEDist>,
                                   Model2GravityEDist>,
                         Model2GravityEDist>;

template class M1SpOptim<System<Model2GravityEDist,
                                ModelEDist>,
                         RankAgent<ToyFeatures5<ModelEDist>,
                                   ModelEDist>,
                         ModelEDist>;

template class M1SpOptim<System<Model2GPowGDist,
                                Model2GPowGDist>,
                         RankAgent<ToyFeatures5<Model2GPowGDist>,
                                   Model2GPowGDist>,
                         Model2GPowGDist>;

template class M1SpOptim<System<Model2EdgeToEdge,
                                Model2EdgeToEdge>,
                         RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                                   Model2EdgeToEdge>,
                         Model2EdgeToEdge>;

template class M1SpOptim<System<Model2EdgeToEdge,
                                ModelIntercept>,
                         RankAgent<ToyFeatures5<ModelIntercept>,
                                   ModelIntercept>,
                         ModelIntercept>;

template class M1SpOptim<System<Model2GPowGDist,
                                ModelGDist>,
                         RankAgent<ToyFeatures5<ModelGDist>,
                                   ModelGDist>,
                         ModelGDist>;

template class M1SpOptim<System<Model2GravityGDist,
                                ModelGDist>,
                         RankAgent<ToyFeatures5<ModelGDist>,
                                   ModelGDist>,
                         ModelGDist>;
