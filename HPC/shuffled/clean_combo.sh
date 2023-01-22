#! /bin/bash
echo "Start!"
Rscript clean_shuf_folder.r -all ./shuffle_30
echo "1/7"
Rscript clean_shuf_folder.r -all ./shuffle_all_30
echo "2/7"
Rscript clean_shuf_folder.r -all ./shuffle_farm_30
echo "3/7"
Rscript clean_shuf_folder.r -all ./shuffle_farm_curveball_30
echo "4/7"
Rscript clean_shuf_folder.r -all ./shuffle_farm_r00_30
echo "5/7"
Rscript clean_shuf_folder.r -all ./shuffle_farm_r0_30
echo "6/7"

Rscript clean_shuf_folder.r -clean ./shuffle_farm_curveball_30_shuff_500
echo "7/7"
echo "Done."
