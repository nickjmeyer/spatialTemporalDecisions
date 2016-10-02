#!/bin/bash
./scripts/run.sh | tee >(sed 's/.*\r//'> run.log)
