#!/bin/bash

cd ../bin

NETS=("rand")
SIZE=("100" "500" "1000")

for j in ${SIZE[@]}
do
    for i in ${NETS[@]}
    do
	./tuneGen "../data/toy/${i}${j}" y;
	submit ./runM1Mles "../data/toy/${i}${j}" y;
	submit ./runM1MlesMiss "../data/toy/${i}${j}" y;
    done
done
