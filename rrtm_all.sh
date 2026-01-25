#! /usr/bin/bash

year=$1
d1=$2
d2=$3

#echo $year

echo "Processing started $(date)"

#mkdir -p /local/storage/RO/rad/pyr/$year

for day in $(seq -f "%03g" $d1 $d2)
do
   python merge_rad.py rrtm $year $day
done

echo "Processing completed $(date)"
