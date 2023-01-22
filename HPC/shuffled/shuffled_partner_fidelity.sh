#! /bin/bash

#$ -q shai.q
#$ -cwd
#$ -N m_shuff_fidelity

##-l h_vmem=2G

LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH

Rscript shuffled_partner_fidelity.R $1