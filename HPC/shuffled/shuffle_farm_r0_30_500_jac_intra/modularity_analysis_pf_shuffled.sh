#! /bin/bash

#$ -q shai.q
#$ -cwd
#$ -N m_shuffled

##-l h_vmem=2G

LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH

Rscript modularity_analysis_pf_shuffled.R $1
