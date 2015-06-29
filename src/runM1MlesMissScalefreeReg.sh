#!/bin/bash

cd ../bin

NETS=("scalefree")
SIZE=("100" "500" "1000")

for i in ${NETS[@]}
do
    for j in ${SIZE[@]}
    do
	submit ./runM1MlesMiss "../data/toy/${i}${j}" y;
    done
done
