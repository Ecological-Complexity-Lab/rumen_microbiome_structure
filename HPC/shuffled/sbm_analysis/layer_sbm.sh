#!/bin/bash

#$ -q shai.q@bhn27
#$ -cwd
#$ -N sbm_layer
#$ -j yes
LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH
Rscript layer_sbm.r $1
