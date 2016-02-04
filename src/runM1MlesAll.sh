#!/bin/bash

cd ../bin

NETS=("grid" "rand" "scalefree")
SIZE=("100" "500" "1000")

for j in ${SIZE[@]}
do
    for i in ${NETS[@]}
    do
	submit ./runM1Mles "../data/toy/${i}${j}" y;
    done
done
