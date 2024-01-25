#!/bin/bash

#$ -q shai.q@bhn27
#$ -cwd
#$ -N embeding
#$ -j yes
LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH
Rscript embed_network.r $1
