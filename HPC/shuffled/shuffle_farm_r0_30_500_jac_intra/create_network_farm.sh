#!/bin/bash

#$ -q shai.q
#$ -cwd
#$ -N Farm
LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH
Rscript create_network_farm.R $1 $2
