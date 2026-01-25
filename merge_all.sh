#! /usr/bin/bash

year=$1
d1=$2
d2=$3

#source /local/storage/miniconda3/etc/profile.d/conda.sh
#conda activate pygeode2

echo $year

for day in $(seq -f "%03g" $d1 $d2)
do
     python merge_rad.py merge $year $day
     python merge_rad.py mergecomp $year $day
done
