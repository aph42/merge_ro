#! /usr/bin/bash

year=$1
m1=$2
m2=$3

#source /local/storage/miniconda3/etc/profile.d/conda.sh
#conda activate pygeode2

echo "Processing started $(date)"

for m in $(seq -f "%02g" $m1 $m2)
do
     #python grid.py sens $year $m
     python grid.py eqrad $year $m
     #python grid.py trop $year $m
done

echo "Processing completed $(date)"
