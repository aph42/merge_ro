
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

# "cosmic2": https://data.cosmic.ucar.edu/gnss-ro/cosmic2/nrt/level2/

declare -a missions=(
            "metopb"
            )

usr='--http-user=phitch --http-passwd=qQSHQ8JhNEjCTY2FUz7x2gHJL'

for mission in "${missions[@]}"
do

for year in {2019..2019}
do

for day in $(seq -f "%03g" 60 365)
do
 
        url='http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data/'$mission'/atmPrf/'$year'.'$day
        fn='/local/storage/RO/dl/'$mission'_atmPrf_'$year'.'$day'.tar'
        path='/local/storage/RO/dl/'$mission'/atmPrf/'$year'.'$day

        echo '==========================='
        echo 'Downloading '$fn
        echo wget -q $usr $url -O $fn
        wget -q $usr $url -O $fn

        echo 'Extracting '$fn 
        tar -xf $fn -C '/local/storage/RO/dl/'

        if [ $? -eq 0 ]
        then
           echo 'Processing '$fn 
           python -m process_mission.py $mission $year $day

           # Clean up downloads
           echo 'Removing downloaded data'
           echo $path
           rm -rf $path
        else
           echo $fn' not valid. Skipping.'
        fi

        echo 'Removing '$fn
        rm -f $fn
done

done

done

