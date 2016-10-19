#!/bin/bash

make -j $(nproc) -f runM1Mles.mk prod
make -j $(nproc) -f runM1MlesMiss.mk prod
make -j $(nproc) -f runM1MlesWns.mk prod
make -j $(nproc) -f runM1MlesWnsMiss.mk prod
make -j $(nproc) -f tuneGen.mk prod
make -j $(nproc) -f tuneGenWNS.mk prod
make -j $(nproc) -f bayesP.mk prod
