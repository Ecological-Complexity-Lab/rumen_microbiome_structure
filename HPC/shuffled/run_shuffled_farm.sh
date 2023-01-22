for i in {1..100}
do
  # assembling needed files in the folder for running network building
  # make folder
  printf -v dir_name "%03d" $i
  rm -r $dir_name 2> /dev/null
  mkdir $dir_name

  #copy sh, r, experiments, functions.r
  cp experiments.csv $dir_name
  cp functions.R $dir_name
  cp create_network_farm.sh $dir_name
  cp create_network_farm.R $dir_name
  
  #  add abundance file
  cp shuff_farm_$dir_name.csv $dir_name
  echo folder $dir_name is ready
  
  #run the jobs in their folder
  cd $dir_name
  qsub create_network_farm.sh $i IT1
  qsub create_network_farm.sh $i IT2
  qsub create_network_farm.sh $i IT3
  qsub create_network_farm.sh $i FI1
  qsub create_network_farm.sh $i UK1
  qsub create_network_farm.sh $i UK2
  qsub create_network_farm.sh $i SE1
  cd ../
done
