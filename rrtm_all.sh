#! /usr/bin/bash

year=$1

echo $year

mkdir -p /local/storage/RO/rad/pyr/$year

for day in $(seq -f "%03g" 1 366)
do
   python merge_rad.py rrtm $year $day
done

