#! /usr/bin/bash

#           "metopa2016"
#           "metopb2016"
#           "sacc"
#           "tsx"
#           "grace"
#           "cnofs"
#           "kompsat5"
#           "cnofsrt"
#           "cosmic2013"
#           "cosmic"
#           "champ2016"
#           "metopb"
#           "metopa"
#           "geoopt"
#           "metopc"
#           "paz"
#           "planetiq"
#           "spire"
#           "tdx"
#           "cosmic2021"
#           "tdx"
#           "geoopt"
#           "metopc"
#           "paz"
#           "planetiq"
#           "spire"

# "cosmic2": https://data.cosmic.ucar.edu/gnss-ro/cosmic2/nrt/level2/
# https://data.cosmic.ucar.edu/gnss-ro/metopa/postProc/level2/2021/331/atmPrf_postProc_2021_331.tar.gz

declare -a missions=(
#           "cosmic2"
#           "tdx"
#           "geoopt"
#           "metopc"
#           "paz"
#           "planetiq"
#           "spire"
            "cosmic2021"
            )

root='/local/storage/RO'

y1=$1
y2=$2

echo $y1 $y2

for mission in "${missions[@]}"
do

for year in $(seq $y1 $y2)
do

for day in $(seq -f "%03g" 1 366)
do
 
        #fn=$root'/dl/'$mission'_atmPrf_'$year'.'$day'.tar'
        mpath=$root'/dl/'$mission
        fn=$mpath'/'$mission'_atmPrf_'$year'.'$day'.tar.gz'
        path=$mpath'/atmPrf/'$year'.'$day

        nc=$root'/raw/'$mission'/'$mission'_atmPrf_'$(date -d "$year-01-01 +$day days - 1 day" "+%Y-%m-%d.nc")
        if [ -f $nc ]
        then
           # Check if modification is 
           modtime=$(stat -c %Y $nc)
           refdate=$(date --date='2026-01-09' +"%s")
           if [ "$modtime" -gt "$refdate" ]
           then
              echo "$nc exists and is more recent than Jan 9th. Skipping."
              continue
           else
              echo "$nc exists but is older than Jan 9th. Re-processing."
           fi
        fi

        echo '==========================='
        echo 'Extracting '$fn 

        mkdir -p $path
        tar -xzf $fn -C $path

        if [ $? -eq 0 ]
        then
           echo 'Processing '$fn 
           python process_mission.py $mission $year $day

           # Clean up unpacked files
           echo 'Removing downloaded data '$path
           rm -rf $path
        else
           echo $fn' not valid. Skipping.'
        fi
done

done

done

