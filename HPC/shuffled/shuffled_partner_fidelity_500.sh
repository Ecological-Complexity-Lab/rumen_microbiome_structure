for j in {1..500}
do
  # move stuff to folder
  printf -v dir_name "%03d" $j
  cp shuffled_partner_fidelity.sh $dir_name
  cp shuffled_partner_fidelity.R $dir_name
  cp fitted_asvs_phylo_tree.rds $dir_name

  cd $dir_name
  qsub shuffled_partner_fidelity.sh $j
  cd ../
done
