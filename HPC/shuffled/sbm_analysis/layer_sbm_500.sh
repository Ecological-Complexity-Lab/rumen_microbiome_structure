#!/bin/bash

# consts:
shuff_folder=../shuffle_farm_r0_30_500_jac_intra
work_folder=../sbm_analysis

# run:
cd $shuff_folder # Go to the shuffled directory to get folder names
args=()
for d in */ 
do
  args+=("$PWD/$d")
done

cd $work_folder # Go to the working directory to run the jobs from
for i in "${args[@]}" # loop on the foldernames to put as job arguments
do
   qsub layer_sbm.sh $i
done