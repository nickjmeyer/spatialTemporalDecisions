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


## Generating network structures

In order to simulate the spread of the disease, a network structure is
needed.  In the `src/` directory is an `R` script to generate
networks.  Generate networks of size `n` use `Rscript genToyNets.R`.
Before generating these networks, make sure a `data` directory exists
in the project root directory.

## Simulating disease with intervention

An example of setting up a run script can be found in
`src/main/runM1Mles.cpp`.  There are multiple models to simulate
disease from.  All model types have their own source and header file.
Model `XYZ`, the files will be named `src/main/modelXYZ.[hc]pp`.  A
model will attempt to load files from a directory inside the generated
network directory.  By default, all parameters are initialized to zero
and a warning is printed if parameter files cannot be found.
