for j in {1..500}
do
  # move stuff to folder
  printf -v dir_name "%03d" $j
  cp Infomap $dir_name
  cp modularity_analysis_shuffled_farm.sh $dir_name
  cp modularity_analysis_shuffled_farm.R $dir_name
  cp core_ASV_30.csv $dir_name # change this according to the observed file used

  cd $dir_name
  #make sure there are execution permission
  chmod +x Infomap
  qsub modularity_analysis_shuffled_farm.sh $j
  cd ../
done

