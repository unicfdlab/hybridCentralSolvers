#!/bin/bash

rm -rf log

blockMesh | tee -a log
pisoCentralFoam | tee -a log
sample | tee -a log
