#!/bin/bash

make -j 48 -f runM1Mles.mk prod
make -j 48 -f runM1MlesMiss.mk prod
make -j 48 -f runM1MlesWns.mk prod
make -j 48 -f runM1MlesWnsMiss.mk prod
make -j 48 -f tuneGen.mk prod
make -j 48 -f tuneGenWNS.mk prod
