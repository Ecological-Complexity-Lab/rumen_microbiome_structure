for j in {1..500}
do
  # move stuff to folder
  printf -v dir_name "%03d" $j
  cp Infomap $dir_name
  cp functions.R $dir_name
  cp modularity_analysis_pf_shuffled.sh $dir_name
  cp modularity_analysis_pf_shuffled.R $dir_name
  cp fitted_asvs_phylo_tree.rds $dir_name
  
  cd $dir_name
  #make sure there are execution permission
  chmod +x Infomap
  qsub modularity_analysis_pf_shuffled.sh $j
  cd ../
done

