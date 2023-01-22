#!/bin/bash

#$ -q shai.q
#$ -cwd
#$ -N Triangles
LD_LIBRARY_PATH=/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/lib/:$LD_LIBRARY_PATH
Rscript farm_triangle_phylogeny.r $1
