for j in {1..50}
do
  # move stuff to folder
  printf -v dir_name "%02d" $j
  mkdir $dir_name
  cp Infomap $dir_name
  cp test_interlayer_amount.R $dir_name
  cp test_interlayer_amount_itr.sh $dir_name
  
  cd $dir_name
  #make sure there are execution permission
  chmod +x Infomap
  qsub test_interlayer_amount_itr.sh $j
  cd ../
done