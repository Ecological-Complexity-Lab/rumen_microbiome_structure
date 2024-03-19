#!/bin/bash

#$ -q shai.q@bhn27
#$ -cwd
#$ -N multiedge
#$ -j yes
LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/lib/:$LD_LIBRARY_PATH
Rscript test_interlayer_amount.R $1
