#!/bin/bash

#$ -q shai.q@bhn27
#$ -cwd
#$ -N Trans
LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH
Rscript Transitivity.R $1
