#include "runner.hpp"



template <class S, class A>
RunStats
TrainRunner<S,A>
::run(S system,
        A agent,
        const int numReps, const int numPoints){
    // double value=0;
    int r,t;
    RunStats rs;
    for(r=0; r<numReps; r++){
        system.model.assignRand(system.paramGen_r,system.paramEst_r);
        system.revert();
        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0)
                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst,system.paramEst);
            system.updateStatus();

            system.nextPoint();
        }
        rs(system.value());
        // value += system.value();
    }
    return rs;
}




template <class S, class A>
RunStats
PlainRunner<S,A>
::run(S system,
        A agent,
        const int numReps, const int numPoints){
    RunStats rs;
    // double value=0;
    int r,t;
    for(r=0; r<numReps; r++){
        // njm::timer.start("plain setup");
        if(system.modelGen_r.sample()){
            std::vector<double> newPar = system.modelGen_r.getPar();
            system.modelEst_r.putPar(newPar.begin());
            system.modelGen_r.setFill(system.sD,system.tD,system.fD,system.dD);
            system.modelEst_r.setFill(system.sD,system.tD,system.fD,system.dD);
        }

        system.revert();
        // njm::timer.stop("plain setup");

        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0){
                // njm::timer.start("plain applyTrt");
                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);
                // njm::timer.stop("plain applyTrt");
            }

            // njm::timer.start("plain updateStatus");
            system.updateStatus();
            // njm::timer.stop("plain updateStatus");

            // njm::timer.start("plain nextPoint");
            system.nextPoint();
            // njm::timer.stop("plain nextPoint");
        }

        rs(system.value());
        // value += system.value();
    }

    return rs;
}

template class PlainRunner<System<ModelGDist,
                                  ModelGDist>,
                           RankAgent<ToyFeatures5<ModelGDist>,
                                     ModelGDist> >;


template class PlainRunner<System<ModelGravityGDist,
                                  ModelGravityGDist>,
                           NoTrt<ModelGravityGDist> >;


template class PlainRunner<System<ModelGravityGDist,
                                  ModelGravityGDist>,
                           RankAgent<WnsFeatures3<ModelGravityGDist>,
                                     ModelGravityGDist> >;


template class PlainRunner<System<Model2GravityGDist,
                                  Model2GravityGDist>,
                           RankAgent<WnsFeatures3<Model2GravityGDist>,
                                     Model2GravityGDist> >;


template class PlainRunner<System<ModelIntercept,
                                  ModelIntercept>,
                           RankAgent<WnsFeatures3<ModelIntercept>,
                                     ModelIntercept> >;


template class PlainRunner<System<Model2EdgeToEdge,
                                  Model2EdgeToEdge>,
                           RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                                     Model2EdgeToEdge> >;


template class PlainRunner<System<Model2GravityEDist,
                                  Model2GravityEDist>,
                           RankAgent<WnsFeatures3<Model2GravityEDist>,
                                     Model2GravityEDist> >;


template class PlainRunner<System<ModelEDist,
                                  ModelEDist>,
                           RankAgent<WnsFeatures3<ModelEDist>,
                                     ModelEDist> >;


template class PlainRunner<System<Model2GravityGDist,
                                  Model2GravityGDist>,
                           RankAgent<ToyFeatures5<Model2GravityGDist>,
                                     Model2GravityGDist> >;


template class PlainRunner<System<Model2GravityEDist,
                                  Model2GravityEDist>,
                           RankAgent<ToyFeatures5<Model2GravityEDist>,
                                     Model2GravityEDist> >;


template class PlainRunner<System<ModelEDist,
                                  ModelEDist>,
                           RankAgent<ToyFeatures5<ModelEDist>,
                                     ModelEDist> >;


template class PlainRunner<System<Model2GPowGDist,
                                  Model2GPowGDist>,
                           RankAgent<ToyFeatures5<Model2GPowGDist>,
                                     Model2GPowGDist> >;


template class PlainRunner<System<Model2EdgeToEdge,
                                  Model2EdgeToEdge>,
                           RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                                     Model2EdgeToEdge> >;


template class PlainRunner<System<ModelIntercept,
                                  ModelIntercept>,
                           RankAgent<ToyFeatures5<ModelIntercept>,
                                     ModelIntercept> >;


template class PlainRunner<System<Model2GPowGDist,
                                  ModelGDist>,
                           RankAgent<ToyFeatures5<ModelGDist>,
                                     ModelGDist> >;

template class PlainRunner<System<ModelGDist,
                                  ModelGDist>,
                           RankAgent<WnsFeatures3<ModelGDist>,
                                     ModelGDist> >;




template<class S, class A>
RunStats
VanillaRunner<S,A>
::run(S system,
        A agent,
        const int numReps, const int numPoints,
        const Starts & starts){
    // double value=0;
    int r,t;

    RunStats rs;
    std::vector<std::vector<double> > valueAll(numReps);
#pragma omp parallel for num_threads(omp_get_max_threads())	\
    shared(valueAll,starts,rs)                              \
    firstprivate(system,agent)                              \
    private(r,t)
    for(r=0; r<numReps; r++){
        njm::resetSeed(r);
        system.reset(starts[r]);

#pragma omp critical
        {
            valueAll.at(r).clear();
            valueAll.at(r).push_back(system.value());
        }
        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0)
                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);

            system.updateStatus();

            system.nextPoint();

#pragma omp critical
            {
                valueAll.at(r).push_back(system.value());
            }
        }

#pragma omp critical
        {
            rs(system.value());
            // value += system.value();
        }

        system.sD.history.push_back(system.sD.status);
        njm::toFile(njm::toString(system.sD.history,"\n","")
                ,njm::sett.datExt(agent.name +
                        "_history_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));
    }
    njm::toFile(njm::toString(valueAll,"\n",""),
            njm::sett.datExt(agent.name+
                    "_values_",".txt"));
    return rs;
}

template class VanillaRunner<System<ModelGravityGDist,
                                    ModelGravityGDist>,
                             NoTrt<ModelGravityGDist> >;

template class VanillaRunner<System<Model2GravityGDist,
                                    Model2GravityGDist>,
                             NoTrt<Model2GravityGDist> >;

template class VanillaRunner<System<Model2GravityEDist,
                                    Model2GravityEDist>,
                             NoTrt<Model2GravityEDist> >;

template class VanillaRunner<System<Model2GPowGDist,
                                    Model2GPowGDist>,
                             NoTrt<Model2GPowGDist> >;

template class VanillaRunner<System<Model2EdgeToEdge,
                                    Model2EdgeToEdge>,
                             NoTrt<Model2EdgeToEdge> >;

template class VanillaRunner<System<Model2GPowGDist,
                                    ModelGDist>,
                             NoTrt<ModelGDist> >;

template class VanillaRunner<System<ModelGravityGDist,
                                    ModelGravityGDist>,
                             ProximalAgent<ModelGravityGDist> >;

template class VanillaRunner<System<Model2GravityGDist,
                                    Model2GravityGDist>,
                             ProximalAgent<Model2GravityGDist> >;


template class VanillaRunner<System<Model2GravityEDist,
                                    Model2GravityEDist>,
                             ProximalAgent<Model2GravityEDist> >;

template class VanillaRunner<System<Model2GPowGDist,
                                    Model2GPowGDist>,
                             ProximalAgent<Model2GPowGDist> >;

template class VanillaRunner<System<Model2EdgeToEdge,
                                    Model2EdgeToEdge>,
                             ProximalAgent<Model2EdgeToEdge> >;

template class VanillaRunner<System<Model2GPowGDist,
                                    ModelGDist>,
                             ProximalAgent<ModelGDist> >;

template class VanillaRunner<System<Model2GravityEDist,
                                    Model2GravityEDist>,
                             AllAgent<Model2GravityEDist> >;

template class VanillaRunner<System<Model2EdgeToEdge,
                                    Model2EdgeToEdge>,
                             AllAgent<Model2EdgeToEdge> >;



template<class S, class A>
RunStats
VanillaRunnerNS<S,A>
::run(S system,
        A agent,
        const int numReps, const int numPoints,
        const Starts & starts){
    // double value=0;
    int r,t;

    std::vector<double> vals(numReps);

    RunStats rs;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
    shared(starts,rs,vals)                                  \
    firstprivate(system,agent)                              \
    private(r,t)
    for(r=0; r<numReps; r++){
        njm::resetSeed(r);
        system.reset(starts[r]);
        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0)
                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);

            system.updateStatus();

            system.nextPoint();

        }

#pragma omp critical
        {
            rs(system.value());
            vals.at(r) = system.value();
            // value += system.value();
        }

    }

    std::cout << std::setprecision(17)
              << std::accumulate(vals.begin(),vals.end(),0.)
              << std::endl;

    RunStats rs2;
    std::for_each(vals.begin(),vals.end(),
            [&rs2](const double & x) {
                rs2(x);
            });
    std::cout << std::setprecision(17)
              << rs2.sMean()
              << std::endl;
    return rs;
}


template class VanillaRunnerNS<System<ModelGravityGDist,
                                      ModelGravityGDist>,
                               NoTrt<ModelGravityGDist> >;

template class VanillaRunnerNS<System<Model2EdgeToEdge,
                                      Model2EdgeToEdge>,
                               NoTrt<Model2EdgeToEdge> >;

template class VanillaRunnerNS<System<Model2GravityEDist,
                                      Model2GravityEDist>,
                               NoTrt<Model2GravityEDist> >;

template class VanillaRunnerNS<System<Model2GravityGDist,
                                      Model2GravityGDist>,
                               NoTrt<Model2GravityGDist> >;

template class VanillaRunnerNS<System<Model2GravityGDist,
                                      Model2GravityGDist>,
                               AllAgent<Model2GravityGDist> >;

template class VanillaRunnerNS<System<Model2EdgeToEdge,
                                      Model2EdgeToEdge>,
                               AllAgent<Model2EdgeToEdge> >;

template class VanillaRunnerNS<System<Model2GravityEDist,
                                      Model2GravityEDist>,
                               AllAgent<Model2GravityEDist> >;

template class VanillaRunnerNS<System<Model2GPowGDist,
                                      Model2GPowGDist>,
                               NoTrt<Model2GPowGDist> >;

template class VanillaRunnerNS<System<Model2GPowGDist,
                                      Model2GPowGDist>,
                               AllAgent<Model2GPowGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
                                      ModelGravityGDist>,
                               ProximalAgent<ModelGravityGDist> >;

template class VanillaRunnerNS<System<Model2GravityGDist,
                                      Model2GravityGDist>,
                               ProximalAgent<Model2GravityGDist> >;

template class VanillaRunnerNS<System<Model2EdgeToEdge,
                                      Model2EdgeToEdge>,
                               ProximalAgent<Model2EdgeToEdge> >;

template class VanillaRunnerNS<System<Model2GravityEDist,
                                      Model2GravityEDist>,
                               ProximalAgent<Model2GravityEDist> >;

template class VanillaRunnerNS<System<Model2GPowGDist,
                                      Model2GPowGDist>,
                               ProximalAgent<Model2GPowGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
                                      ModelGravityGDist>,
                               MyopicAgent<ModelGravityGDist> >;

template class VanillaRunnerNS<System<Model2GravityGDist,
                                      Model2GravityGDist>,
                               MyopicAgent<Model2GravityGDist> >;

template class VanillaRunnerNS<System<Model2EdgeToEdge,
                                      Model2EdgeToEdge>,
                               MyopicAgent<Model2EdgeToEdge> >;

template class VanillaRunnerNS<System<Model2GravityEDist,
                                      Model2GravityEDist>,
                               MyopicAgent<Model2GravityEDist> >;

template class VanillaRunnerNS<System<Model2GPowGDist,
                                      Model2GPowGDist>,
                               MyopicAgent<Model2GPowGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
                                      ModelGravityGDist>,
                               RankAgent<ToyFeatures5<ModelGravityGDist>,
                                         ModelGravityGDist> >;

template class VanillaRunnerNS<System<ModelGravityGDist,
                                      ModelGravityGDist>,
                               RankAgent<WnsFeatures3<ModelGravityGDist>,
                                         ModelGravityGDist> >;

template class VanillaRunnerNS<System<Model2GravityGDist,
                                      Model2GravityGDist>,
                               RankAgent<WnsFeatures3<Model2GravityGDist>,
                                         Model2GravityGDist> >;

template class VanillaRunnerNS<System<Model2EdgeToEdge,
                                      Model2EdgeToEdge>,
                               RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                                         Model2EdgeToEdge> >;

template class VanillaRunnerNS<System<Model2GravityEDist,
                                      Model2GravityEDist>,
                               RankAgent<WnsFeatures3<Model2GravityEDist>,
                                         Model2GravityEDist> >;

template class VanillaRunnerNS<System<Model2GravityGDist,
                                      Model2GravityGDist>,
                               RankAgent<ToyFeatures5<Model2GravityGDist>,
                                         Model2GravityGDist> >;

template class VanillaRunnerNS<System<Model2EdgeToEdge,
                                      Model2EdgeToEdge>,
                               RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                                         Model2EdgeToEdge> >;

template class VanillaRunnerNS<System<Model2GravityEDist,
                                      Model2GravityEDist>,
                               RankAgent<ToyFeatures5<Model2GravityEDist>,
                                         Model2GravityEDist> >;

template class VanillaRunnerNS<System<Model2GPowGDist,
                                      Model2GPowGDist>,
                               RankAgent<ToyFeatures5<Model2GPowGDist>,
                                         Model2GPowGDist> >;



template <class S, class A>
RunStats
FitOnlyRunner<S,A>
::run(S system,
        A agent,
        const int numReps, const int numPoints,
        const Starts & starts){
    // double value=0;
    int r,t;

    RunStats rs;
    std::vector<std::vector<double> > valueAll(numReps);
    std::vector< std::vector<double> > pars;
#pragma omp parallel for num_threads(omp_get_max_threads())	\
    shared(valueAll,starts,rs)                      \
    firstprivate(pars,system,agent)                 \
    private(r,t)
    for(r=0; r<numReps; r++){
        njm::resetSeed(r);
        system.reset(starts[r]);

#pragma omp critical
        {
            valueAll.at(r).clear();
            valueAll.at(r).push_back(system.value());
        }
        pars.clear();

        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0){
                system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
                        t > system.fD.trtStart);

                // only save pars after estimating
                pars.push_back(system.modelEst.getPar());

                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);
            }

            system.updateStatus();

            system.nextPoint();

#pragma omp critical
            {
                valueAll.at(r).push_back(system.value());
            }
        }

#pragma omp critical
        {
            rs(system.value());
            // value += system.value();
        }

        system.sD.history.push_back(system.sD.status);
        njm::toFile(njm::toString(system.sD.history,"\n","")
                ,njm::sett.datExt(agent.name+
                        "_history_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));

        // write estimated pars to file
        njm::toFile(njm::toString(pars,"\n","")
                ,njm::sett.datExt(agent.name+
                        "_pars_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));

    }
    njm::toFile(njm::toString(valueAll,"\n",""),
            njm::sett.datExt(agent.name+
                    "_values_",".txt"));
    return rs;
}


template class FitOnlyRunner<System<ModelGravityGDist,
                                    ModelGravityGDist>,
                             MyopicAgent<ModelGravityGDist> >;

template class FitOnlyRunner<System<Model2GravityGDist,
                                    Model2GravityGDist>,
                             MyopicAgent<Model2GravityGDist> >;

template class FitOnlyRunner<System<Model2GravityEDist,
                                    Model2GravityEDist>,
                             MyopicAgent<Model2GravityEDist> >;

template class FitOnlyRunner<System<Model2GravityEDist,
                                    ModelEDist>,
                             MyopicAgent<ModelEDist> >;

template class FitOnlyRunner<System<Model2GPowGDist,
                                    Model2GPowGDist>,
                             MyopicAgent<Model2GPowGDist> >;

template class FitOnlyRunner<System<Model2EdgeToEdge,
                                    Model2EdgeToEdge>,
                             MyopicAgent<Model2EdgeToEdge> >;

template class FitOnlyRunner<System<Model2EdgeToEdge,
                                    ModelIntercept>,
                             MyopicAgent<ModelIntercept> >;

template class FitOnlyRunner<System<Model2GPowGDist,
                                    ModelGDist>,
                             MyopicAgent<ModelGDist> >;

template class FitOnlyRunner<System<Model2GravityGDist,
                                    ModelGDist>,
                             MyopicAgent<ModelGDist> >;



template <class S, class A, class Optim>
RunStats
OptimRunner<S,A,Optim>
::run(S system,
        A agent,
        Optim optim,
        const int numReps, const int numPoints,
        const Starts & starts){
    int tick,tickR,tock,tockR,done=0;
    tick = std::time(NULL);
    double hours;

    RunStats rs;
    // double value=0;
    int r,t;
    std::vector<std::vector<double> > valueAll(numReps);
    std::vector<std::vector<double> > pars;
    std::vector<std::vector<double> > weights;

    std::vector<int> times(numReps);

    // int threads = (omp_get_max_threads() < 16 ? 1 : omp_get_max_threads());
    int threads = omp_get_max_threads();

#pragma omp parallel for num_threads(threads)     \
    shared(valueAll,tock,tick,starts,rs)          \
    firstprivate(system,agent,optim,pars,weights) \
    private(r,t,tockR,tickR)
    for(r=0; r<numReps; r++){
        njm::resetSeed(r);
        // record time for each replication
        tickR=std::time(NULL);

        system.reset(starts[r]);
        agent.reset();
        optim.reset();

#pragma omp critical
        {
            valueAll.at(r).clear();
            valueAll.at(r).push_back(system.value());
        }
        pars.clear();
        weights.clear();

        // begin rep r
        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0){
                // njm::timer.start("fit");
                system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
                        t > system.fD.trtStart);
                // njm::timer.stop("fit");

                // only save pars after estimating
                pars.push_back(system.modelEst.getPar());

                // njm::timer.start("optim");
                optim.optim(system,agent);
                // njm::timer.stop("optim");

                weights.push_back(agent.tp.getPar());

                // njm::timer.start("trt");
                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);
                // njm::timer.stop("trt");
            }

            system.updateStatus();

            system.nextPoint();

#pragma omp critical
            {
                valueAll.at(r).push_back(system.value());
            }
        }
        // end rep r...time to write results to disk

#pragma omp critical
        {
            rs(system.value());
            // value += system.value();
        }

        // store time for iteration
        tockR = std::time(NULL);
#pragma omp critical
        {
            times.at(r) = tockR - tickR;
        }

        // write history to file
        system.sD.history.push_back(system.sD.status);
        njm::toFile(njm::toString(system.sD.history,"\n","")
                ,njm::sett.datExt(agent.name+"_"+optim.name+
                        "_history_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));

        // write estimated pars to file
        njm::toFile(njm::toString(pars,"\n","")
                ,njm::sett.datExt(agent.name+"_"+optim.name+
                        "_pars_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));

        // write weights
        njm::toFile(njm::toString(weights,"\n","")
                ,njm::sett.datExt(agent.name+"_"+optim.name+
                        "_weights_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));

        // write optim parameters to file
        njm::toFile(njm::toString(optim.tp.getPar()," ","\n")
                ,njm::sett.datExt(agent.name+"_"+optim.name+
                        "_tunePar_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));


#pragma omp critical
        {
            done++;
            tock = std::time(NULL);
            hours = ((double)(tock-tick))/((double)done);
            hours /= 3600.0;

            njm::toFile("Completed " + njm::toString(done,"",6,0) +
                    " out of " + njm::toString(numReps,"",6,0) +
                    " in " + njm::toString(hours,"",8,4) + " hours" +
                    " with value " + njm::toString(rs.sMean(),"",6,4) +
                    "\n",
                    njm::sett.datExt(agent.name+"_"+optim.name+"_status_",
                            ".txt"));
        }

    }

    njm::toFile(njm::toString(valueAll,"\n",""),
            njm::sett.datExt(agent.name+"_"+optim.name+
                    "_values_",".txt"));

    njm::toFile(njm::toString(times,"\n",""),
            njm::sett.datExt(agent.name+"_"+optim.name+
                    "_times_",".txt"));


    return rs;
}






template class
OptimRunner<System<Model2GravityGDist,
                   Model2GravityGDist>,
            RankAgent<ToyFeatures5<Model2GravityGDist>,
                      Model2GravityGDist>,
            M1SpOptim<System<Model2GravityGDist,
                             Model2GravityGDist>,
                      RankAgent<ToyFeatures5<Model2GravityGDist>,
                                Model2GravityGDist>,
                      Model2GravityGDist> >;

template class
OptimRunner<System<Model2GravityEDist,
                   Model2GravityEDist>,
            RankAgent<ToyFeatures5<Model2GravityEDist>,
                      Model2GravityEDist>,
            M1SpOptim<System<Model2GravityEDist,
                             Model2GravityEDist>,
                      RankAgent<ToyFeatures5<Model2GravityEDist>,
                                Model2GravityEDist>,
                      Model2GravityEDist> >;


template class
OptimRunner<System<Model2GravityEDist,
                   ModelEDist>,
            RankAgent<ToyFeatures5<ModelEDist>,
                      ModelEDist>,
            M1SpOptim<System<Model2GravityEDist,
                             ModelEDist>,
                      RankAgent<ToyFeatures5<ModelEDist>,
                                ModelEDist>,
                      ModelEDist> >;


template class
OptimRunner<System<Model2GPowGDist,
                   Model2GPowGDist>,
            RankAgent<ToyFeatures5<Model2GPowGDist>,
                      Model2GPowGDist>,
            M1SpOptim<System<Model2GPowGDist,
                             Model2GPowGDist>,
                      RankAgent<ToyFeatures5<Model2GPowGDist>,
                                Model2GPowGDist>,
                      Model2GPowGDist> >;


template class
OptimRunner<System<Model2EdgeToEdge,
                   Model2EdgeToEdge>,
            RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                      Model2EdgeToEdge>,
            M1SpOptim<System<Model2EdgeToEdge,
                             Model2EdgeToEdge>,
                      RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                                Model2EdgeToEdge>,
                      Model2EdgeToEdge> >;


template class
OptimRunner<System<Model2EdgeToEdge,
                   ModelIntercept>,
            RankAgent<ToyFeatures5<ModelIntercept>,
                      ModelIntercept>,
            M1SpOptim<System<Model2EdgeToEdge,
                             ModelIntercept>,
                      RankAgent<ToyFeatures5<ModelIntercept>,
                                ModelIntercept>,
                      ModelIntercept> >;


template class
OptimRunner<System<Model2GPowGDist,
                   ModelGDist>,
            RankAgent<ToyFeatures5<ModelGDist>,
                      ModelGDist>,
            M1SpOptim<System<Model2GPowGDist,
                             ModelGDist>,
                      RankAgent<ToyFeatures5<ModelGDist>,
                                ModelGDist>,
                      ModelGDist> >;


template class
OptimRunner<System<Model2GravityGDist,
                   ModelGDist>,
            RankAgent<ToyFeatures5<ModelGDist>,
                      ModelGDist>,
            M1SpOptim<System<Model2GravityGDist,
                             ModelGDist>,
                      RankAgent<ToyFeatures5<ModelGDist>,
                                ModelGDist>,
                      ModelGDist> >;


template class
OptimRunner<System<Model2EdgeToEdge,
                   ModelIntercept>,
            RankAgent<WnsFeatures3<ModelIntercept>,
                      ModelIntercept>,
            M1SpOptim<System<Model2EdgeToEdge,
                             ModelIntercept>,
                      RankAgent<WnsFeatures3<ModelIntercept>,
                                ModelIntercept>,
                      ModelIntercept> >;


template class
OptimRunner<System<Model2GravityEDist,
                   ModelEDist>,
            RankAgent<WnsFeatures3<ModelEDist>,
                      ModelEDist>,
            M1SpOptim<System<Model2GravityEDist,
                             ModelEDist>,
                      RankAgent<WnsFeatures3<ModelEDist>,
                                ModelEDist>,
                      ModelEDist> >;


template class
OptimRunner<System<Model2EdgeToEdge,
                   Model2EdgeToEdge>,
            RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                      Model2EdgeToEdge>,
            M1SpOptim<System<Model2EdgeToEdge,
                             Model2EdgeToEdge>,
                      RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                                Model2EdgeToEdge>,
                      Model2EdgeToEdge> >;


template class
OptimRunner<System<Model2GravityGDist,
                   Model2GravityGDist>,
            RankAgent<WnsFeatures3<Model2GravityGDist>,
                      Model2GravityGDist>,
            M1SpOptim<System<Model2GravityGDist,
                             Model2GravityGDist>,
                      RankAgent<WnsFeatures3<Model2GravityGDist>,
                                Model2GravityGDist>,
                      Model2GravityGDist> >;

template class
OptimRunner<System<Model2GravityEDist,
                   Model2GravityEDist>,
            RankAgent<WnsFeatures3<Model2GravityEDist>,
                      Model2GravityEDist>,
            M1SpOptim<System<Model2GravityEDist,
                             Model2GravityEDist>,
                      RankAgent<WnsFeatures3<Model2GravityEDist>,
                                Model2GravityEDist>,
                      Model2GravityEDist> >;

template class
OptimRunner<System<Model2GravityGDist,
                   ModelGDist>,
            RankAgent<WnsFeatures3<ModelGDist>,
                      ModelGDist>,
            M1SpOptim<System<Model2GravityGDist,
                             ModelGDist>,
                      RankAgent<WnsFeatures3<ModelGDist>,
                                ModelGDist>,
                      ModelGDist> >;



template <class S, class A, class Optim>
RunStats
OptimRunnerNS<S,A,Optim>
::run(S system,
        A agent,
        Optim optim,
        const int numReps, const int numPoints,
        const Starts & starts){

    RunStats rs;
    // double value=0;
    int r,t;

    int threads = omp_get_max_threads();

#pragma omp parallel for num_threads(threads)   \
    shared(starts,rs)                           \
    firstprivate(system,agent,optim)            \
    private(r,t)
    for(r=0; r<numReps; r++){
        system.reset(starts[r]);
        agent.reset();
        optim.reset();

        // begin rep r
        for(t=system.sD.time; t<numPoints; t++){
            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0){
                system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
                        t > system.fD.trtStart);

                optim.optim(system,agent);

                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);
            }

            system.updateStatus();

            system.nextPoint();

        }
        // end rep r

#pragma omp critical
        {
            rs(system.value());
            // value += system.value();
        }

    }

    return rs;
}



template <class S, class A, class Optim>
RunStats
TuneRunner<S,A,Optim>
::run(S system,
        A agent,
        Optim optim,
        const int numReps, const int numPoints){

    RunStats rs;
    // double value=0;
    int r,t;
    for(r=0; r<numReps; r++){
        if(system.modelGen_r.sample()){
            std::vector<double> newPar = system.modelGen_r.getPar();
            system.modelEst_r.putPar(newPar.begin());
        }

        system.revert();
        agent.reset();

        // begin rep r
        for(t=system.sD.time; t<numPoints; t++){

            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0){
                system.modelEst.fit(system.sD,system.tD,system.fD,system.dD,
                        t > system.fD.trtStart);

                optim.optim(system,agent);

                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst);
            }

            system.updateStatus();

            system.nextPoint();

        }
        // end rep r

        rs(system.value());
        // value += system.value();
    }

    return rs;
}


template class
TuneRunner<System<ModelGDist,
                  ModelGDist>,
           RankAgent<ToyFeatures5<ModelGDist>,
                     ModelGDist>,
           M1SpOptim<System<ModelGDist,
                            ModelGDist>,
                     RankAgent<ToyFeatures5<ModelGDist>,
                               ModelGDist>,
                     ModelGDist> >;


template class
TuneRunner<System<ModelGDist,
                  ModelGDist>,
           RankAgent<WnsFeatures3<ModelGDist>,
                     ModelGDist>,
           M1SpOptim<System<ModelGDist,
                            ModelGDist>,
                     RankAgent<WnsFeatures3<ModelGDist>,
                               ModelGDist>,
                     ModelGDist> >;


template class
TuneRunner<System<ModelGravityGDist,
                  ModelGravityGDist>,
           RankAgent<WnsFeatures3<ModelGravityGDist>,
                     ModelGravityGDist>,
           M1SpOptim<System<ModelGravityGDist,
                            ModelGravityGDist>,
                     RankAgent<WnsFeatures3<ModelGravityGDist>,
                               ModelGravityGDist>,
                     ModelGravityGDist> >;


template class
TuneRunner<System<Model2GravityGDist,
                  Model2GravityGDist>,
           RankAgent<WnsFeatures3<Model2GravityGDist>,
                     Model2GravityGDist>,
           M1SpOptim<System<Model2GravityGDist,
                            Model2GravityGDist>,
                     RankAgent<WnsFeatures3<Model2GravityGDist>,
                               Model2GravityGDist>,
                     Model2GravityGDist> >;

template class
TuneRunner<System<ModelIntercept,
                  ModelIntercept>,
           RankAgent<WnsFeatures3<ModelIntercept>,
                     ModelIntercept>,
           M1SpOptim<System<ModelIntercept,
                            ModelIntercept>,
                     RankAgent<WnsFeatures3<ModelIntercept>,
                               ModelIntercept>,
                     ModelIntercept> >;


template class
TuneRunner<System<Model2EdgeToEdge,
                  Model2EdgeToEdge>,
           RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                     Model2EdgeToEdge>,
           M1SpOptim<System<Model2EdgeToEdge,
                            Model2EdgeToEdge>,
                     RankAgent<WnsFeatures3<Model2EdgeToEdge>,
                               Model2EdgeToEdge>,
                     Model2EdgeToEdge> >;


template class
TuneRunner<System<Model2GravityEDist,
                  Model2GravityEDist>,
           RankAgent<WnsFeatures3<Model2GravityEDist>,
                     Model2GravityEDist>,
           M1SpOptim<System<Model2GravityEDist,
                            Model2GravityEDist>,
                     RankAgent<WnsFeatures3<Model2GravityEDist>,
                               Model2GravityEDist>,
                     Model2GravityEDist> >;


template class
TuneRunner<System<ModelEDist,
                  ModelEDist>,
           RankAgent<WnsFeatures3<ModelEDist>,
                     ModelEDist>,
           M1SpOptim<System<ModelEDist,
                            ModelEDist>,
                     RankAgent<WnsFeatures3<ModelEDist>,
                               ModelEDist>,
                     ModelEDist> >;


template class
TuneRunner<System<Model2GravityGDist,
                  Model2GravityGDist>,
           RankAgent<ToyFeatures5<Model2GravityGDist>,
                     Model2GravityGDist>,
           M1SpOptim<System<Model2GravityGDist,
                            Model2GravityGDist>,
                     RankAgent<ToyFeatures5<Model2GravityGDist>,
                               Model2GravityGDist>,
                     Model2GravityGDist> >;


template class
TuneRunner<System<Model2GravityEDist,
                  Model2GravityEDist>,
           RankAgent<ToyFeatures5<Model2GravityEDist>,
                     Model2GravityEDist>,
           M1SpOptim<System<Model2GravityEDist,
                            Model2GravityEDist>,
                     RankAgent<ToyFeatures5<Model2GravityEDist>,
                               Model2GravityEDist>,
                     Model2GravityEDist> >;


template class
TuneRunner<System<ModelEDist,
                  ModelEDist>,
           RankAgent<ToyFeatures5<ModelEDist>,
                     ModelEDist>,
           M1SpOptim<System<ModelEDist,
                            ModelEDist>,
                     RankAgent<ToyFeatures5<ModelEDist>,
                               ModelEDist>,
                     ModelEDist> >;


template class
TuneRunner<System<Model2GPowGDist,
                  Model2GPowGDist>,
           RankAgent<ToyFeatures5<Model2GPowGDist>,
                     Model2GPowGDist>,
           M1SpOptim<System<Model2GPowGDist,
                            Model2GPowGDist>,
                     RankAgent<ToyFeatures5<Model2GPowGDist>,
                               Model2GPowGDist>,
                     Model2GPowGDist> >;


template class
TuneRunner<System<Model2EdgeToEdge,
                  Model2EdgeToEdge>,
           RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                     Model2EdgeToEdge>,
           M1SpOptim<System<Model2EdgeToEdge,
                            Model2EdgeToEdge>,
                     RankAgent<ToyFeatures5<Model2EdgeToEdge>,
                               Model2EdgeToEdge>,
                     Model2EdgeToEdge> >;


template class
TuneRunner<System<ModelIntercept,
                  ModelIntercept>,
           RankAgent<ToyFeatures5<ModelIntercept>,
                     ModelIntercept>,
           M1SpOptim<System<ModelIntercept,
                            ModelIntercept>,
                     RankAgent<ToyFeatures5<ModelIntercept>,
                               ModelIntercept>,
                     ModelIntercept> >;


template class
TuneRunner<System<Model2GPowGDist,
                  ModelGDist>,
           RankAgent<ToyFeatures5<ModelGDist>,
                     ModelGDist>,
           M1SpOptim<System<Model2GPowGDist,
                            ModelGDist>,
                     RankAgent<ToyFeatures5<ModelGDist>,
                               ModelGDist>,
                     ModelGDist> >;




template <class S, class A, class Optim>
RunStats
TestRunner<S,A,Optim>
::run(S system,
        A agent,
        Optim optim,
        const int numReps, const int numPoints){

    RunStats rs;
    // double value=0;
    int r,t;
    std::vector<std::vector<double> > valueAll(numReps);
    std::vector<std::vector<double> > weights;
    int threads = (omp_get_max_threads() < 16 ? 1 : omp_get_max_threads());
#pragma omp parallel for num_threads(threads)   \
    shared(valueAll,rs)                         \
    firstprivate(system,agent,optim,weights)    \
    private(r,t)
    for(r=0; r<numReps; r++){
        system.reset();
#pragma omp critical
        {
            valueAll.at(r).clear();
            valueAll.at(r).push_back(system.value());
        }
        weights.clear();
        for(t=system.sD.time; t<numPoints; t++){
            if(t==system.fD.trtStart && system.sD.numNotInfec > 0){
                optim.optim(system,agent);
                weights.push_back(agent.tp.getPar());
            }

            if(t>=system.fD.trtStart && system.sD.numNotInfec > 0)
                agent.applyTrt(system.sD,system.tD,system.fD,system.dD,
                        system.modelEst,system.paramEst);

            system.updateStatus();

            system.nextPoint();

#pragma omp critical
            {
                valueAll.at(r).push_back(system.value());
            }
        }

#pragma omp critical
        {
            rs(system.value());
            // value += system.value();
        }

        system.sD.history.push_back(system.sD.status);
        njm::toFile(njm::toString(system.sD.history,"\n","")
                ,njm::sett.datExt(agent.name+"_"+optim.name+
                        "_history_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));
        njm::toFile(njm::toString(weights,"\n","")
                ,njm::sett.datExt(agent.name+"_"+optim.name+
                        "_weights_"+
                        njm::toString(r,"",0,0)
                        +"_",".txt"));
    }
    njm::toFile(njm::toString(valueAll,"\n",""),
            njm::sett.datExt(agent.name+"_"+optim.name+
                    "_values_",".txt"));
    return rs;
}
