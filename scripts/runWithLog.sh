#!/bin/bash
./scripts/run.sh 2>&1 | tee >(sed 's/.*\r//'> run.log)
