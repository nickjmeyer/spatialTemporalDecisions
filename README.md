# Optimal treatment allocations in space and time for online control of emerging infectious diseases

Simulate the spread of a disease under multiple intervention strategies.

## Compiling the code

To compile the simulation code, run the following.
- `mkdir build && cd build`
- `cmake ..`
- `make`

To run test code execute `make test`.

Libraries needed to compile the code are
- GSL
- Armadillo
- Boost (system, filesystem, thread)
- Google Logging
- Google Flags
- Google Test

## Setting up White Nose Syndrome (WNS) data and toy network structures

Place all WNS data in `data/wns`.  All toy network structures will be created by scripts and placed in `data/toy/<name of network structure>`.

### WNS data

The data files are
- `fips.txt`: Fips code by county.
- `eDist.txt`: Euclidean distance.
- `gDist.txt`: Geodesic distance.
- `caves.txt`: Number of caves per county.
- `xcov.txt`: Covariates for each county (details are in the paper).
- `network.txt`: Adjacency matrix for counties.
- `centroidsLong.txt`: Longitude of county centroid.
- `centroidsLat.txt`: Latitude of county centroid.
- `subGraph.txt`: Subgraph centrality.
- `betweenness.txt`: Betweenness centrality.
- `trtStart.txt`: Time point to start treatment.
- `startingLocations.txt`: Index of initial infected locations.
- `period.txt`: A legacy value used to reduce computational cost of the stochastic optimization algorithm.  This value governed the simulation horizon when estimating the objective function.
- `obsData.txt`: Observed WNS infections.  Each row is a time point and each column is a county.  For any value `x`, if `floor(x / 2) == 1` the location is infected.  The value `2` indicates infection because in the simulations `x % 2 == 1` indicates treatment.

### Generating network structures

In order to simulate the spread of the disease, a network structure is
needed.  In the `src/` directory is an `R` script to generate
networks.  Generate networks of size `n` run `Rscript genToyNets.R n`.
Before generating these networks, make sure a `data` directory exists
in the project root directory.

## Fitting models

Each model can be fit using the executables that start with `mcmc`.  This will fit to the WNS data.  To convert the model to the toy networks, the parameters have to be copied to the toy network directory (use the `copyParams` executable) then tuned to the network (use the `tuneGen` executable).  Examples of this process can be found in the `scripts/individualRunners`.

## Simulating disease with intervention

An example of setting up a run script can be found in
`src/main/runM1Mles.cpp`.  There are multiple models to simulate
disease from.  All model types have their own source and header file.
Model `XYZ`, the files will be named `src/main/modelXYZ.[hc]pp`.  A
model will attempt to load files from a directory inside the generated
network directory.  By default, all parameters are initialized to zero
and a warning is printed if parameter files cannot be found.

## Running simulations

The results from the paper were run using the following scripts
- `runCrpSpatialCorrLog.sh`
- `runCrpSpatialMissLog.sh`
- `runGridEdgeCorrLog.sh`
- `runGridEdgeMissLog.sh`
- `runGridSpatialCorrLog.sh`
- `runGridSpatialMissLog.sh`
- `runRandEdgeCorrLog.sh`
- `runRandEdgeMissLog.sh`
- `runRandSpatialCorrLog.sh`
- `runRandSpatialMissLog.sh`
- `runScalefreeEdgeCorrLog.sh`
- `runScalefreeEdgeMissLog.sh`
- `runWnsEdgeCorrLog.sh`
- `runWnsEdgeMissLog.sh`
- `runWnsSpatialCorrLog.sh`
- `runWnsSpatialMissLog.sh`
