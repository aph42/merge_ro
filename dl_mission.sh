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
#           "cosmic2021"
#           "cosmic2"
#           "tdx"
#           "geoopt"
#           "metopc"
#           "paz"
#           "planetiq"
            "spire"
            )

usr='--http-user=phitch --http-passwd=qQSHQ8JhNEjCTY2FUz7x2gHJL'
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
        # Old URL and filename w/ authentication
        #url='http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data/'$mission'/atmPrf/'$year'.'$day
        #fn=$root'/dl/'$mission'_atmPrf_'$year'.'$day'.tar'

        # New path for many missions
        #url='https://data.cosmic.ucar.edu/gnss-ro/'$mission'/postProc/level2/'$year'/'$day'/atmPrf_postProc_'$year'_'$day'.tar.gz'

        # geoopt: post processed (to 2021-259?; use NRT for now)
        #url='https://data.cosmic.ucar.edu/gnss-ro/geoopt/noaa/postProc/level2/'$year'/'$day'/atmPrf_postProc_'$year'_'$day'.tar.gz'

        # planetiq/geoopt/spire: near-real-time (all dates?)
        url='https://data.cosmic.ucar.edu/gnss-ro/'$mission'/noaa/nrt/level2/'$year'/'$day'/atmPrf_nrt_'$year'_'$day'.tar.gz'

        # Cosmic 2 path
        #url='https://data.cosmic.ucar.edu/gnss-ro/'$mission'/nrt/level2/'$year'/'$day'/atmPrf_nrt_'$year'_'$day'.tar.gz'

        # Cosmic 1 2021 reprocessing
        #url='https://data.cosmic.ucar.edu/gnss-ro/cosmic1/repro2021/level2/'$year'/'$day'/atmPrf_repro2021_'$year'_'$day'.tar.gz'

        mpath=$root'/dl/'$mission
        fn=$mpath'/'$mission'_atmPrf_'$year'.'$day'.tar.gz'
        path=$mpath'/'$mission'/atmPrf/'$year'.'$day


        if [ -f $fn ]
        then
           echo "$fn exists. Skipping."
           continue
        fi

        mkdir -p $mpath

        echo '==========================='
        echo 'Downloading '$fn
        echo wget $url -O $fn
        #echo wget -q $usr $url -O $fn
        #wget -q $usr $url -O $fn
        wget $url -O $fn

        if [ $? -eq 0 ]
        then
           echo 'Downloaded '$fn 
        else
           echo $fn' not valid.'
           rm -f $fn
        fi

        #echo 'Removing '$fn
        #rm -f $fn
done

done

done

