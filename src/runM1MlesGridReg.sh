#!/bin/bash

cd ../bin

NETS=("grid")
SIZE=("100" "500" "1000" "10000")

for i in ${NETS[@]}
do
    for j in ${SIZE[@]}
    do
	./tuneGen "../data/toy/${i}${j}" y;
	submit ./runM1Mles "../data/toy/${i}${j}" y;
	submit ./runM1MlesMiss "../data/toy/${i}${j}" y;
    done
done
